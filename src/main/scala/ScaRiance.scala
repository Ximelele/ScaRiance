import org.apache.spark.SparkContext
import org.apache.spark.sql.SparkSession

import scala.collection.parallel.CollectionConverters.*

case class ScaRiance(control_file: String, tumour_file: String, skip_allele_counting: Boolean = false, skip_imputation: Boolean = false, skip_segmentation: Boolean = false):

  private val utils = Utils()
  private val prepare_Wgs = PrepareWgs()
  private val impute = Impute()
  private val haplotype = Haplotype()
  private val spark: SparkSession = SparkSession.builder()
    .appName("ScaRiance")
    .master("local[*]")
    .getOrCreate()

  def run(): Unit = {


    this.setDefaultValues(spark)

    if (!skip_allele_counting) {
      alle_counting()
    }

    if (!skip_imputation) {
      imputation()
    }

    if (!skip_segmentation) {
      segmentation()
    }

    spark.stop()
  }


  def alle_counting(): Unit = {
    prepare_Wgs.prepareWgs(spark = spark, utils = utils, controlFile = control_file, tumourFile = tumour_file)
  }

  def imputation(): Unit = {
    this.utils.chromosomeNames.foreach(chrom => {

      impute.runHaplotyping(spark, chrom.toString, this.utils)
      haplotype.getChromosomeBafs(spark = spark, SNP_file = s"${utils.allele_directory}/${utils.tumourName}_alleleFrequencies_$chrom.txt", haplotype_File = s"${utils.impute_directory}/${utils.tumourName}_impute_output_${chrom}_allHaplotypeInfo.txt", utils = utils, output_file = s"${utils.impute_directory}/${utils.tumourName}_impute_output_${chrom}_heterozygousMutBAFs_haplotyped.txt", minCounts = 10)
    })
    utils.concatenateBAFfiles(spark = spark, inputStart = s"${utils.impute_directory}/${utils.tumourName}_impute_output_")
  }

  def segmentation(): Unit = {

    def runRcode(cmd: Seq[String]): Unit = {
      import scala.sys.process.*
      val exitCode = cmd.!

      // Check the exit code
      if (exitCode == 0) {
        println("R script executed successfully")
      } else {
        println(s"R script failed with exit code $exitCode")
      }
    }

    val cmd_segmentation = Seq(
      "Rscript",
      "-e",
      s"""
              source("/app/ScalaBattenberg/src/main/R/Segmentation.R")

              segment.baf.phased(
                samplename = "${utils.tumourName}",
                inputfile = "${utils.impute_directory}/${utils.tumourName}_heterozygousMutBAFs_haplotyped.txt",
                outputfile = "${utils.working_directory}/${utils.tumourName}.BAFsegmented.txt",
              )
              """
    )

    runRcode(cmd_segmentation)


    // fit copy number
    val logr_file = s"${utils.working_directory}/${utils.tumourName}_mutantLogR.tab"

    val input_baf_segment = s"${utils.working_directory}/${utils.tumourName}.BAFsegmented.txt"
    val input_baf = s"${utils.working_directory}/${utils.tumourName}_mutantBAF.tab"

    val outputfile_prefix = s"${utils.working_directory}/${utils.tumourName}_"
    val ploting_prefix = s"${utils.plots_directory}/${utils.tumourName}_"
    val log_segment_file = s"${utils.impute_directory}/${utils.tumourName}.logRsegmented.txt"
    val cmd_fit_copy_number = Seq(
      "Rscript",
      "-e",
      s"""
                  source("/app/ScalaBattenberg/src/main/R/fitCopyNumber.R")

                  fit.copy.number(
                    samplename = "${utils.tumourName}",
                    outputfile.prefix  = "$outputfile_prefix",
                    inputfile.baf.segmented = "$input_baf_segment",
                    inputfile.baf = "$input_baf",
                    inputfile.logr = "$logr_file",
                    log_segment_file = "$log_segment_file",
                    ploting_prefix = "$ploting_prefix"
                  )
                  """
    )

    runRcode(cmd_fit_copy_number)

    val call_subclones = Seq(
      "Rscript",
      "-e",
      s"""
                  source("/app/ScalaBattenberg/src/main/R/fitCopyNumber.R")

                  callSubclones(
                    sample.name = "${utils.tumourName}",
                    baf.segmented.file = "$input_baf_segment",
                    logr.file = "$logr_file",
                    rho.psi.file = "${utils.working_directory}/${utils.tumourName}_rho_and_psi.txt",
                    output.file = "${utils.working_directory}/${utils.tumourName}_copynumber.txt",
                    output.gw.figures.prefix= "${utils.plots_directory}/${utils.tumourName}_BattenbergProfile",
                    masking_output_file="${utils.plots_directory}/${utils.tumourName}_segment_masking_details.txt",
                    chr_names = c(${utils.chromosomeNames.map(_.toString).mkString("\"", "\", \"", "\"")})

                  )
                  """
    )

    runRcode(call_subclones)


    val callChrXsubclones = Seq(
      "Rscript",
      "-e",
      s"""
                      source("/app/ScalaBattenberg/src/main/R/fitCopyNumber.R")

                      callChrXsubclones(
                        tumourname = "${utils.working_directory}/${utils.tumourName}",
                        chrom_names = c(${utils.chromosomeNames.map(_.toString).mkString("\"", "\", \"", "\"")})

                      )
                      """
    )

    runRcode(callChrXsubclones)
  }


  private def setDefaultValues(spark: SparkSession): Unit = {
    this.utils.tumourName = tumour_file.split("/").last
    this.utils.controlName = control_file.split("/").last
    this.utils.is_male = this.utils.isMale(tumour_file, spark)
    this.utils.setChromosomesNames(control_file)

    if (this.utils.chromosomeNames(0).toString.contains("chr")) {
      this.utils.referenciesFile.setPostfixDirectory()
    } else {
      this.utils.referenciesFile.setPostFixFile()
    }

    this.utils.createPatientDirectory(this.utils.controlName)

  }


