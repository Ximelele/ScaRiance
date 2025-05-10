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

    // allele skipping is happening inside function
    alle_counting()

    if (!skip_imputation) {
      imputation()
    }

    if (!skip_segmentation) {
      segmentation()
    }

    val ploting_prefix = s"${utils.plots_directory}/${utils.tumourName}_"
    this.utils.chromosomeNames.foreach(chrom => {
      val plots = Seq(
        "Rscript",
        "-e",
        s"""
                              source("/app/ScalaBattenberg/src/main/R/plotting.R")

                              plot.haplotype.data(
                                haplotyped.baf.file = "${utils.impute_directory}/${utils.tumourName}_chr${chrom.stripPrefix("chr")}_heterozygousMutBAFs_haplotyped.txt",
                                imageFileName = "${ploting_prefix}_chr${chrom.stripPrefix("chr")}_heterozygousData.png",
                                chrom = "$chrom")
                              )
                              """
      )
    })

    spark.stop()
  }


  def alle_counting(): Unit = {
    prepare_Wgs.prepareWgs(spark = spark, utils = utils, controlFile = control_file, tumourFile = tumour_file, skip_alle_counting = skip_allele_counting)
  }

  def imputation(): Unit = {
    this.utils.chromosomeNames.foreach(chrom => {

      impute.runHaplotyping(spark, chrom.toString, this.utils)
      haplotype.getChromosomeBafs(spark = spark, SNP_file = s"${utils.allele_directory}/${utils.tumourName}_alleleFrequencies_$chrom.txt", haplotype_File = s"${utils.impute_directory}/${utils.tumourName}_impute_output_${chrom}_allHaplotypeInfo.txt", utils = utils, output_file = s"${utils.impute_directory}/${utils.tumourName}_impute_output_${chrom}_heterozygousMutBAFs_haplotyped.txt", minCounts = 10)
    })
    utils.concatenateBAFfiles(spark = spark, inputStart = s"${utils.impute_directory}/${utils.tumourName}_impute_output_")
  }

  def runRcode(cmd: Seq[String]): Unit = {
    import scala.sys.process.*
    val exitCode = cmd.!

    // Check the exit code
    if (exitCode == 0) {
      println("R script executed successfully")
    } else {
      println(s"R script failed with exit code $exitCode")
      System.exit(1)
    }
  }

  def segmentation(): Unit = {



    // add print segmenting
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
    // add print fit copy
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

    // add sublonal
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


