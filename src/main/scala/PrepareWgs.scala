import org.apache.spark.sql.SparkSession

import scala.sys.process.*

case class PrepareWgs():

  def prepareWgs(utils: Utils, controlFile: String, tumourFile: String): Unit = {

    var tumourAlleleCountsFilePrefix: String = s"${utils.allele_directory}/${utils.tumourName}_alleleFrequencies_"
    var normalAlleleCountsFilePrefix: String = s"${utils.allele_directory}/${utils.controlName}_alleleFrequencies_"
    //    utils.chromosomeNames.foreach(chromsome => {
    //      getAlleleCounts(bamFile = tumourFile, outputFile = s"${utils.allele_directory}/${tumourFile.split("/").last}_alleleFrequencies_$chromsome.txt", g1000Loci = s"${utils.g1000prefix}$chromsome.txt")
    //
    //      getAlleleCounts(bamFile = controlFile, outputFile = s"${utils.allele_directory}/${controlFile.split("/").last}_alleleFrequencies_$chromsome.txt", g1000Loci = s"${utils.g1000prefix}$chromsome.txt")
    //    })

    val spark = SparkSession.builder()
      .appName("Battenberg")
      .master("local[*]")
      .getOrCreate()

    utils.processAlleleData(spark = spark, tumourAlleleCountsFilePrefix = tumourAlleleCountsFilePrefix, normalAlleleCountsFilePrefix = normalAlleleCountsFilePrefix, minCounts = 10, samplename = utils.tumourName, BAFnormalFile = s"${utils.working_directory}/${utils.tumourName}_normalBAF.tab", BAFmutantFile = s"${utils.working_directory}/${utils.tumourName}_mutantBAF.tab", logRnormalFile = s"${utils.working_directory}/${utils.tumourName}_normalLogR.tab", logRmutantFile = s"${utils.working_directory}/${utils.tumourName}_mutantLogR.tab", combinedAlleleCountsFile = s"${utils.working_directory}/${utils.tumourName}_alleleCounts.tab")

  }

  def getAlleleCounts(bamFile: String, outputFile: String, g1000Loci: String, minBaseQual: Int = 20, minMapQual: Int = 35, alleleCounter: String = "alleleCounter"): Unit = {
    //    var alleleCounterCommand = Seq(
    //      alleleCounter,
    //      "-b", bamFile,
    //      "-l", g1000Loci,
    //      "-o", outputFile,
    //      "-m", minBaseQual.toString,
    //      "-q", minMapQual.toString
    //    ).mkString(" ")
    //
    //    val counterVersion: String = s"$alleleCounter --version".!!.trim
    //
    //    if counterVersion.substring(0, 1).toInt >= 4 then
    //      alleleCounterCommand = s"$alleleCounterCommand --dense-snps"
    //
    //    val exitCode: Int = alleleCounterCommand.!
    //
    //    require(exitCode == 0, s"Command failed with exit code $exitCode")
  }