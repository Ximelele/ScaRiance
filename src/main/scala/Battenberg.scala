import org.apache.spark.sql.SparkSession

import scala.collection.parallel.CollectionConverters.*

case class Battenberg(control_file: String, tumour_file: String):

  private val utils = Utils()
  private val prepare_Wgs = PrepareWgs()
  private val impute = Impute()
  private val haplotype = Haplotype()

  def run(): Unit = {
    this.setDefaultValues()


    //    this.male = this.utils.isMale(tumour_file)
    val spark = SparkSession.builder()
      .appName("Battenberg")
      .master("local[*]")
      .getOrCreate()
    //    prepare_Wgs.prepareWgs(utils = utils, controlFile = control_file, tumourFile = tumour_file)
    //    println("Macka Macka")
    //
    //    this.utils.chromosomeNames.foreach(chrom => {
    //      println(chrom)
    //      //impute.runHaplotyping(spark, chrom, this.utils)
    //      haplotype.getChromosomeBafs(spark = spark, SNP_file = s"${utils.allele_directory}/${utils.tumourName}_alleleFrequencies_$chrom.txt", haplotype_File = s"${utils.impute_directory}/${utils.tumourName}_impute_output_${chrom}_allHaplotypeInfo.txt", utils = utils, output_file = s"${utils.impute_directory}/${utils.tumourName}_impute_output_${chrom}_heterozygousMutBAFs_haplotyped.txt", minCounts = 10)
    //    })
    //    utils.concatenateBAFfiles(spark = spark, inputStart = s"${utils.impute_directory}/${utils.tumourName}_impute_output_")


    // Starting R segmetation part
    import scala.sys.process._
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

    var exitCode = cmd_segmentation.!

    // Check the exit code
    if (exitCode == 0) {
      println("R script executed successfully")
    } else {
      println(s"R script failed with exit code $exitCode")
    }

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

    exitCode = cmd_fit_copy_number.!

    // Check the exit code
    if (exitCode == 0) {
      println("R script executed successfully")
    } else {
      println(s"R script failed with exit code $exitCode")
    }

    val call_sublocnes = Seq(
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
                output.gw.figures.prefix=paste(samplename[sampleidx],"_BattenbergProfile", sep=""),
                masking_output_file=paste(samplename[sampleidx], "_segment_masking_details.txt", sep=""),
                chr_names= c(${utils.chromosomeNames.mkString(",")})

              )
              """
    )

    exitCode = call_sublocnes.!

    if (exitCode == 0) {
      println("R script executed successfully")
    } else {
      println(s"R script failed with exit code $exitCode")
    }

  }

  private def setDefaultValues(): Unit = {
    this.utils.tumourName = tumour_file.split("/").last
    this.utils.controlName = control_file.split("/").last
    this.utils.is_male = this.utils.isMale(tumour_file)
    this.utils.setChromosomesNames(control_file)

    if (this.utils.chromosomeNames(0).toString.contains("chr")) {
      this.utils.referenciesFile.setPostfixDirectory()
    } else {
      this.utils.referenciesFile.setPostFixFile()
    }

    this.utils.createPatientDirectory(this.utils.controlName)

  }


