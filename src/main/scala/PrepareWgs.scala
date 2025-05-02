import org.apache.spark.sql.*
import org.apache.spark.sql.functions.*
import org.apache.spark.sql.types.*
import org.apache.spark.storage.StorageLevel

import java.io.File
import java.nio.file.{Files, Paths, StandardCopyOption}
import scala.collection.mutable.ListBuffer
import scala.collection.parallel.CollectionConverters.*
import scala.collection.parallel.ParSeq
import scala.sys.process.*
import scala.util.{Random, Try}

case class PrepareWgs():

  def prepareWgs(spark: SparkSession, utils: Utils, controlFile: String, tumourFile: String): Unit = {

    val tumourAlleleCountsFilePrefix: String = s"${utils.allele_directory}/${utils.tumourName}_alleleFrequencies_"
    val normalAlleleCountsFilePrefix: String = s"${utils.allele_directory}/${utils.controlName}_alleleFrequencies_"
    utils.chromosomeNames.foreach(chromosome => {
      alleleCounting(bamFile = tumourFile, outputFile = s"$tumourAlleleCountsFilePrefix$chromosome.txt", g1000Loci = s"${utils.referenciesFile.g1000prefix}$chromosome.txt")

      alleleCounting(bamFile = controlFile, outputFile = s"$normalAlleleCountsFilePrefix$chromosome.txt", g1000Loci = s"${utils.referenciesFile.g1000prefix}$chromosome.txt")
    })


    processAlleleData(spark = spark, tumourAlleleCountsFilePrefix = tumourAlleleCountsFilePrefix, normalAlleleCountsFilePrefix = normalAlleleCountsFilePrefix, minCounts = 10, samplename = utils.tumourName, BAFnormalFile = s"${utils.working_directory}/${utils.tumourName}_normalBAF.tab", BAFmutantFile = s"${utils.working_directory}/${utils.tumourName}_mutantBAF.tab", logRnormalFile = s"${utils.working_directory}/${utils.tumourName}_normalLogR.tab", logRmutantFile = s"${utils.working_directory}/${utils.tumourName}_mutantLogR.tab", combinedAlleleCountsFile = s"${utils.working_directory}/${utils.tumourName}_alleleCounts.tab", utils = utils)

  }

  private def alleleCounting(bamFile: String, outputFile: String, g1000Loci: String, minBaseQual: Int = 20, minMapQual: Int = 35, alleleCounter: String = "alleleCounter"): Unit = {
    var alleleCounterCommand = Seq(
      alleleCounter,
      "-b", bamFile,
      "-l", g1000Loci,
      "-o", outputFile,
      "-m", minBaseQual.toString,
      "-q", minMapQual.toString
    ).mkString(" ")

    val counterVersion: String = s"$alleleCounter --version".!!.trim

    if counterVersion.substring(0, 1).toInt >= 4 then
      alleleCounterCommand = s"$alleleCounterCommand --dense-snps"

    val exitCode: Int = alleleCounterCommand.!

    require(exitCode == 0, s"Command failed with exit code $exitCode")
  }

  private def processAlleleData(
                                 spark: SparkSession,
                                 tumourAlleleCountsFilePrefix: String,
                                 normalAlleleCountsFilePrefix: String,
                                 minCounts: Int = 10,
                                 samplename: String,
                                 BAFnormalFile: String,
                                 BAFmutantFile: String,
                                 logRnormalFile: String,
                                 logRmutantFile: String,
                                 combinedAlleleCountsFile: String,
                                 utils: Utils
                               ): Unit = {


    var inputData = utils.concatenateAlleleCountFiles(spark, tumourAlleleCountsFilePrefix).persist(StorageLevel.MEMORY_AND_DISK)
    var normalInputData = utils.concatenateAlleleCountFiles(spark, normalAlleleCountsFilePrefix).persist(StorageLevel.MEMORY_AND_DISK)
    var alleleData = utils.concatenateG1000SnpFiles(spark, utils.referenciesFile.g1000alleleprefix).persist(StorageLevel.MEMORY_AND_DISK)


    alleleData = alleleData
      .withColumnRenamed("CHR", "allele_CHR")
      .withColumnRenamed("position", "allele_POS")
      .withColumn("chrpos_allele", concat(col("allele_CHR"), lit("_"), col("allele_POS")))

    normalInputData = normalInputData
      .withColumnRenamed("#CHR", "normal_CHR")
      .withColumnRenamed("POS", "normal_POS")
      .withColumnRenamed("Count_A", "normal_Count_A")
      .withColumnRenamed("Count_C", "normal_Count_C")
      .withColumnRenamed("Count_G", "normal_Count_G")
      .withColumnRenamed("Count_T", "normal_Count_T")
      .withColumn("chrpos_normal", concat(col("normal_CHR"), lit("_"), col("normal_POS")))
      .drop("Good_depth")

    inputData = inputData
      .withColumnRenamed("#CHR", "mutant_CHR")
      .withColumnRenamed("POS", "mutant_POS")
      .withColumnRenamed("Count_A", "mutant_Count_A")
      .withColumnRenamed("Count_C", "mutant_Count_C")
      .withColumnRenamed("Count_G", "mutant_Count_G")
      .withColumnRenamed("Count_T", "mutant_Count_T")
      .withColumn("chrpos_tumour", concat(col("mutant_CHR"), lit("_"), col("mutant_POS")))
      .drop("Good_depth")

    val extractChr = regexp_extract(col("chrpos"), "^(?:chr)?(\\d+|X)_.*", 1)
    val extractPos = regexp_extract(col("chrpos"), "^(?:chr)?(?:\\d+|X)_(\\d+)", 1).cast("int")

    var joinedDf = alleleData.join(normalInputData, alleleData("chrpos_allele") === normalInputData("chrpos_normal"),
        "inner")
      .join(inputData, alleleData("chrpos_allele") === inputData("chrpos_tumour"),
        "inner")

    inputData.unpersist()
    inputData = null

    alleleData.unpersist()
    alleleData = null

    normalInputData.unpersist()
    normalInputData = null


    joinedDf = joinedDf.drop("chrpos_allele")
      .drop("chrpos_normal")
      .drop("allele_POS")
      .drop("allele_CHR")
      .drop("normal_CHR")
      .drop("normal_POS")
      .withColumnRenamed("chrpos_tumour", "chrpos")
      .orderBy(
        when(extractChr === "X", 100).otherwise(extractChr.cast("int")),
        extractPos
      )


    val selectAlleleCount = udf((alleleIdx: Int, countA: Int, countC: Int, countG: Int, countT: Int) => {
      alleleIdx match {
        case 1 => countA
        case 2 => countC
        case 3 => countG
        case 4 => countT
      }
    }, IntegerType)

    joinedDf = joinedDf
      .withColumn("normCount1",
        selectAlleleCount(
          col("a0"),
          col("normal_Count_A"),
          col("normal_Count_C"),
          col("normal_Count_G"),
          col("normal_Count_T")
        )
      )
      .withColumn("normCount2",
        selectAlleleCount(
          col("a1"),
          col("normal_Count_A"),
          col("normal_Count_C"),
          col("normal_Count_G"),
          col("normal_Count_T")
        )
      )
      .withColumn("totalNormal", col("normCount1") + col("normCount2"))
      .withColumn("mutCount1",
        selectAlleleCount(
          col("a0"),
          col("mutant_Count_A"),
          col("mutant_Count_C"),
          col("mutant_Count_G"),
          col("mutant_Count_T")
        )
      )
      .withColumn("mutCount2",
        selectAlleleCount(
          col("a1"),
          col("mutant_Count_A"),
          col("mutant_Count_C"),
          col("mutant_Count_G"),
          col("mutant_Count_T")
        )
      )
      .withColumn("totalmutant", col("mutCount1") + col("mutCount2"))


    joinedDf = joinedDf.drop("chrpos")
      .drop("normal_Count_A", "normal_Count_C", "normal_Count_G", "normal_Count_T")
      .drop("mutant_Count_A", "mutant_Count_C", "mutant_Count_G", "mutant_Count_T")
      .drop("chrpos")
      .drop("a0", "a1")


    joinedDf = if (minCounts > 0) {
      println(s"minCount=$minCounts")
      joinedDf.filter(
        col("totalNormal") >= minCounts && col("totalMutant") >= 1)
    } else {
      joinedDf
    }


    joinedDf = joinedDf.withColumn("selector", (rand() * 2).cast("int"))

    //    val log2UDF = udf((value: Double) => {
    //      scala.math.log(value) / scala.math.log(2.0)
    //    }, DoubleType)


    joinedDf = joinedDf.withColumn("normalBAF",
        when(col("selector") === 0, col("normCount1") / col("totalNormal"))
          .otherwise(col("normCount2") / col("totalNormal"))
      )

      .withColumn("mutantBAF",
        when(col("selector") === 0, col("mutCount1") / col("totalMutant"))
          .otherwise(col("mutCount2") / col("totalMutant"))
      )

      .withColumn("normalLogR", lit(0).cast("int"))
      .withColumn("mutantLogR", col("totalMutant") / col("totalNormal"))

    val meanMutantLogR = joinedDf.agg(avg("mutantLogR")).first().getDouble(0)

    joinedDf = joinedDf.withColumn("log2mutantLogR", log2(col("mutantLogR") / lit(meanMutantLogR)))

    joinedDf = joinedDf.withColumnRenamed("mutant_CHR", "Chromosome")
      .withColumnRenamed("mutant_POS", "Position")


    joinedDf = joinedDf.persist(StorageLevel.DISK_ONLY)

    utils.saveSingleFile(joinedDf.select(col("Chromosome"), col("Position"), col("normalBAF").as(utils.tumourName)), BAFnormalFile)
    utils.saveSingleFile(joinedDf.select(col("Chromosome"), col("Position"), col("mutantBAF").as(utils.tumourName)), BAFmutantFile)
    utils.saveSingleFile(joinedDf.select(col("Chromosome"), col("Position"), col("normalLogR").as(utils.tumourName)), logRnormalFile)
    utils.saveSingleFile(joinedDf.select(col("Chromosome"), col("Position"), col("log2mutantLogR").as(utils.tumourName)), logRmutantFile)
    utils.saveSingleFile(joinedDf.select(col("Chromosome"), col("Position"), col("mutCount1").as("mutCountT1"), col("mutCount2").as("mutCountT2"), col("normCount1").as("mutCountN1"), col("normCount2").as("mutCountN2")), combinedAlleleCountsFile)

    joinedDf.unpersist()
  }