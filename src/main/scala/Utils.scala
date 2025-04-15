import org.apache.spark.sql.{Column, SaveMode}

import java.io.File
import java.nio.file.{Files, Paths, StandardCopyOption}
import scala.collection.mutable.ListBuffer
import scala.collection.parallel.CollectionConverters.*
import scala.collection.parallel.ParSeq
import scala.sys.process._
import scala.util.Try

case class Utils():
  val g1000prefix: String = "/app/references38/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_"
  val g1000alleleprefix: String = "/app/references38/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_allele_index_"
  var working_directory: String = sys.props("user.dir")
  var impute_directory: String = "Impute"
  var plots_directory: String = "Plots"
  var allele_directory: String = "Allele_frequencies"
  var chromosomeNames: ParSeq[Any] = ParSeq()
  var tumourName: String = ""
  var controlName: String = ""

  import org.apache.spark.sql.functions.*
  import org.apache.spark.sql.types.*
  import org.apache.spark.sql.{DataFrame, SparkSession, Row}

  import scala.util.Random

  private def saveSingleFile(df: DataFrame, outputPath: String): Unit = {

    val tempDir = s"${outputPath}_temp"
    df.coalesce(1).write
      .option("header", "true")
      .option("delimiter", "\t")
      .mode(SaveMode.Overwrite)
      .format("csv")
      .save(tempDir)

    val tempDirFile = new File(tempDir)
    val partFile = tempDirFile.listFiles()
      .find(file => file.getName.startsWith("part-") && !file.getName.endsWith(".crc"))
      .getOrElse(throw new RuntimeException(s"No part file found in $tempDir"))

    // Create the output file's parent directory if it doesn't exist
    val outputFile = new File(outputPath)
    Option(outputFile.getParentFile).foreach(_.mkdirs())

    // Copy the part file to the output path
    Files.copy(
      partFile.toPath,
      Paths.get(outputPath),
      StandardCopyOption.REPLACE_EXISTING
    )

    // Delete the temporary directory recursively
    def deleteDirectory(file: File): Unit = {
      if (file.isDirectory) {
        file.listFiles.foreach(deleteDirectory)
      }
      file.delete()
    }

    deleteDirectory(tempDirFile)
  }

  def processAlleleData(
                         spark: SparkSession,
                         tumourAlleleCountsFilePrefix: String,
                         normalAlleleCountsFilePrefix: String,
                         minCounts: Int = 10,
                         samplename: String,
                         BAFnormalFile: String,
                         BAFmutantFile: String,
                         logRnormalFile: String,
                         logRmutantFile: String,
                         combinedAlleleCountsFile: String
                       ): Unit = {

    val seed = System.currentTimeMillis()

    import org.apache.spark.storage.StorageLevel
    val random = new Random(seed)


    var inputData = concatenateAlleleCountFiles(spark, tumourAlleleCountsFilePrefix)
    var normalInputData = concatenateAlleleCountFiles(spark, normalAlleleCountsFilePrefix)
    var alleleData = concatenateG1000SnpFiles(spark, g1000alleleprefix)


    // Create chromosome_position identifiers for joining

    var alleleDataForJoin = alleleData
      .withColumnRenamed("CHR", "allele_CHR")
      .withColumnRenamed("position", "allele_POS")
    alleleData.unpersist()
    alleleData = null
    var normalWithChr = normalInputData
      .withColumnRenamed("#CHR", "normal_CHR")
      .withColumnRenamed("POS", "normal_POS")
    normalInputData.unpersist()
    normalInputData = null

    var inputWithChr = inputData
      .withColumnRenamed("#CHR", "mutant_CHR")
      .withColumnRenamed("POS", "mutant_POS")
    inputData.unpersist()
    inputData = null

    // Join the datasets on CHR and POS
    var alleleDataWithChrPos = alleleDataForJoin
      .withColumn("chrpos_allele", concat(col("allele_CHR"), lit("_"), col("allele_POS")))
    alleleDataForJoin.unpersist()
    alleleDataForJoin = null
    var normalWithChrPos = normalWithChr
      .withColumn("chrpos_normal", concat(col("normal_CHR"), lit("_"), col("normal_POS")))
    normalWithChr.unpersist()
    normalWithChr = null
    var inputWithChrPos = inputWithChr
      .withColumn("chrpos_tumour", concat(col("mutant_CHR"), lit("_"), col("mutant_POS")))
    inputWithChr.unpersist()
    inputWithChr = null
    // Collect the distinct chrpos values from each dataframe
    val alleleChromsPos = alleleDataWithChrPos.select("chrpos_allele").distinct()
    val normalChromsPos = normalWithChrPos.select("chrpos_normal").distinct()
    val tumourChromsPos = inputWithChrPos.select("chrpos_tumour").distinct()

    // Find the intersection
    val extractChr = regexp_extract(col("chrpos"), "^(?:chr)?(\\d+|X)_.*", 1)
    val extractPos = regexp_extract(col("chrpos"), "^(?:chr)?(?:\\d+|X)_(\\d+)", 1).cast("int")


    var intersection1 = alleleChromsPos.join(normalChromsPos,
        alleleChromsPos("chrpos_allele") === normalChromsPos("chrpos_normal"),
        "inner")
      .select("chrpos_allele")
      .withColumnRenamed("chrpos_allele", "chrpos")
      .orderBy(
        // Handle X chromosome specially (X comes after numbers)
        when(extractChr === "X", 100).otherwise(extractChr.cast("int")),
        extractPos
      )


    // Filter each dataframe to only include the matched chromosome-position pairs
    val filteredAlleleData = alleleDataWithChrPos
      .join(intersection1, alleleDataWithChrPos("chrpos_allele") === intersection1("chrpos"), "inner")
      .orderBy(
        // Handle X chromosome specially (X comes after numbers)
        when(extractChr === "X", 100).otherwise(extractChr.cast("int")),
        extractPos
      )
      .drop("chrpos", "chrpos_allele")
    alleleDataWithChrPos.unpersist()
    alleleDataWithChrPos = null

    val filteredNormalData = normalWithChrPos
      .join(intersection1, normalWithChrPos("chrpos_normal") === intersection1("chrpos"), "inner")
      .orderBy(
        // Handle X chromosome specially (X comes after numbers)
        when(extractChr === "X", 100).otherwise(extractChr.cast("int")),
        extractPos
      )
      .drop("chrpos", "chrpos_normal")
    normalWithChrPos.unpersist()
    normalWithChrPos = null

    val filteredInputData = inputWithChrPos
      .join(intersection1, inputWithChrPos("chrpos_tumour") === intersection1("chrpos"), "inner")
      .orderBy(
        // Handle X chromosome specially (X comes after numbers)
        when(extractChr === "X", 100).otherwise(extractChr.cast("int")),
        extractPos
      )
      .drop("chrpos", "chrpos_tumour")

    inputWithChrPos.unpersist()
    inputWithChrPos = null
    intersection1.unpersist()
    intersection1 = null

    val normalData = filteredNormalData.select(
      filteredNormalData.columns.slice(2, 6).map(col): _*
    )
    val mutantData = filteredInputData.select(
      filteredInputData.columns.slice(2, 6).map(col): _*
    )


    import org.apache.spark.sql.expressions.Window
    import org.apache.spark.sql.functions.*
    import org.apache.spark.sql.types.IntegerType
    import spark.implicits.*

    // Extract both allele columns separately
    val a0Indices = filteredAlleleData.select("a0").as[Int].collect()
    val a1Indices = filteredAlleleData.select("a1").as[Int].collect()

    // Create a single UDF that takes the allele type as a parameter
    val selectAllele = udf((rowIdx: Int, alleleType: String) => {
      if (alleleType == "a0") a0Indices(rowIdx) else a1Indices(rowIdx)
    }, IntegerType)

    // Add row index to normal data
    val normalDataWithIndex = normalData.withColumn("row_idx",
      row_number().over(Window.orderBy(monotonically_increasing_id())) - 1)
    val mutantWithIndex = normalData.withColumn("row_idx",
      row_number().over(Window.orderBy(monotonically_increasing_id())) - 1)

    // Create a function for the count selection logic
    def createCountColumn(data: DataFrame, alleleType: String, outputColName: String) = {
      data.withColumn(
        outputColName,
        when(selectAllele(col("row_idx"), lit(alleleType)) === 1, col("Count_A"))
          .when(selectAllele(col("row_idx"), lit(alleleType)) === 2, col("Count_C"))
          .when(selectAllele(col("row_idx"), lit(alleleType)) === 3, col("Count_G"))
          .when(selectAllele(col("row_idx"), lit(alleleType)) === 4, col("Count_T"))
      ).select(outputColName)
    }

    // Use the function to create both columns
    val normCount1 = createCountColumn(normalDataWithIndex, "a0", "normCount1")
    val normCount2 = createCountColumn(normalDataWithIndex, "a1", "normCount2")
    val totalNormalDF = normCount1.join(normCount2).withColumn("total", col("normCount1") + col("normCount2"))
      .select("total")
    val mutCount1 = createCountColumn(mutantWithIndex, "a0", "mutCount1")
    val mutCount2 = createCountColumn(mutantWithIndex, "a1", "mutCount2")
    val totalMutantDF = mutCount1.join(mutCount2).withColumn("total", col("mutCount1") + col("mutCount2"))
      .select("total")

    saveSingleFile(totalNormalDF, "totalNormal")
    saveSingleFile(totalMutantDF, "totalMutant")
    println("Si v pici brasko")

    System.exit(0)


    import org.apache.spark.sql.functions.{lit, when}


    // Extract allele-specific counts and calculate BAF and LogR

    // todo issue bude pravdepodobne tu ako sa handluje getAlleleCount
    //    val processedData = joinedData
    //      // Add random selector for allele swapping
    //      .withColumn("selector", when(rand(seed) > 0.5, 1).otherwise(0))
    //
    //      // Extract counts for normal data
    //      .withColumn("normCount1", getAlleleCount(
    //        col("a0"),
    //        col("normal_Count_A"), col("normal_Count_C"),
    //        col("normal_Count_G"), col("normal_Count_T")
    //      ))
    //      .withColumn("normCount2", getAlleleCount(
    //        col("a1"),
    //        col("normal_Count_A"), col("normal_Count_C"),
    //        col("normal_Count_G"), col("normal_Count_T")
    //      ))
    //
    //      // Extract counts for mutant data (from inputData)
    //      .withColumn("mutCount1", getAlleleCount(
    //        col("a0"),
    //        col("mutant_Count_A"), col("mutant_Count_C"),
    //        col("mutant_Count_G"), col("mutant_Count_T")
    //      ))
    //      .withColumn("mutCount2", getAlleleCount(
    //        col("a1"),
    //        col("mutant_Count_A"), col("mutant_Count_C"),
    //        col("mutant_Count_G"), col("mutant_Count_T")
    //      ))
    //
    //      // Calculate totals
    //      .withColumn("totalNormal", col("normCount1") + col("normCount2"))
    //      .withColumn("totalMutant", col("mutCount1") + col("mutCount2"))
    System.exit(0)
    //    val filteredByCount = if (minCounts > 0) {
    //      println(s"minCount=$minCounts")
    //      processedData.filter(col("totalNormal") >= minCounts && col("totalMutant") >= 1)
    //    } else {
    //      processedData
    //    }
    //
    //    // Calculate BAF and LogR
    //    val finalData = filteredByCount
    //      // Calculate BAF based on selector
    //      .withColumn("normalBAF",
    //        when(col("selector") === 0, col("normCount1") / col("totalNormal"))
    //          .otherwise(col("normCount2") / col("totalNormal"))
    //      )
    //      .withColumn("mutantBAF",
    //        when(col("selector") === 0, col("mutCount1") / col("totalMutant"))
    //          .otherwise(col("mutCount2") / col("totalMutant"))
    //      )
    //      // Calculate LogR values
    //      .withColumn("normalLogR", lit(0))
    //      .withColumn("rawMutantLogR", col("totalMutant") / col("totalNormal"))
    //
    //    // Calculate mean of mutantLogR for normalization
    //    val mutantLogRMean = finalData.agg(mean("rawMutantLogR")).first().getDouble(0)
    //
    //    // Create final datasets for output
    //
    //    val germlineBAF = finalData.select(
    //      col("allele_CHR").as("Chromosome"),
    //      col("allele_POS").as("Position"),
    //      col("normalBAF").as(tumourName)
    //    )
    //
    //    val germlineLogR = finalData.select(
    //      col("allele_CHR").as("Chromosome"),
    //      col("allele_POS").as("Position"),
    //      col("normalLogR").as(tumourName)
    //    )
    //
    //    val tumorBAF = finalData.select(
    //      col("allele_CHR").as("Chromosome"),
    //      col("allele_POS").as("Position"),
    //      col("mutantBAF").as(tumourName)
    //    )
    //
    //    val tumorLogR = finalData.select(
    //      col("allele_CHR").as("Chromosome"),
    //      col("allele_POS").as("Position"),
    //      log2(col("rawMutantLogR") / lit(mutantLogRMean)).as(tumourName)
    //    )
    //
    //    val alleleCounts = finalData.select(
    //      col("allele_CHR").as("Chromosome"),
    //      col("allele_POS").as("Position"),
    //      col("mutCount1"), col("mutCount2"), col("normCount1"), col("normCount2")
    //    )
    //
    //
    //    saveSingleFile(germlineBAF, BAFnormalFile)
    //    // Save data to disk
    //    //    germlineBAF.coalesce(1).write
    //    //      .option("header", "true")
    //    //      .option("delimiter", "\t")
    //    //      .mode("overwrite")
    //    //      .format("csv")
    //    //      .save(BAFnormalFile)
    //
    //    //    tumorBAF.coalesce(1).write
    //    //      .option("header", "true")
    //    //      .option("delimiter", "\t")
    //    //      .mode("overwrite")
    //    //      .format("csv")
    //    //      .save(BAFmutantFile)
    //    saveSingleFile(tumorBAF, BAFmutantFile)
    //
    //    //    germlineLogR.coalesce(1).write
    //    //      .option("header", "true")
    //    //      .option("delimiter", "\t")
    //    //      .mode("overwrite")
    //    //      .format("csv")
    //    //      .save(logRnormalFile)
    //    saveSingleFile(germlineLogR, logRnormalFile)
    //
    //    //    tumorLogR.coalesce(1).write
    //    //      .option("header", "true")
    //    //      .option("delimiter", "\t")
    //    //      .mode("overwrite")
    //    //      .format("csv")
    //    //      .save(logRmutantFile)
    //
    //    saveSingleFile(tumorLogR, logRmutantFile)
    //
    //    //    alleleCounts.coalesce(1).write
    //    //      .option("header", "true")
    //    //      .option("delimiter", "\t")
    //    //      .mode("overwrite").format("csv").save(combinedAlleleCountsFile)
    //
    //    saveSingleFile(alleleCounts, combinedAlleleCountsFile)


  }

  /**
   * Generic reading function tailored for reading in genomic data using Spark
   *
   * @param spark    The SparkSession instance
   * @param file     Filename of the file to read in
   * @param header   Whether the file contains a header (Default: true)
   * @param rowNames Whether the file contains row names (Default: false)
   * @param sep      Column separator (Default: "\t")
   * @param chromCol The column number that contains chromosome denominations. This column will automatically be cast as a String. Should be counted including the rowNames (Default: 1)
   * @param skip     The number of rows to skip before reading (Default: 0)
   * @return A DataFrame with contents of the file
   */
  def readTableGeneric(
                        spark: SparkSession,
                        file: String,
                        header: Boolean = true,
                        rowNames: Boolean = false,
                        sep: String = "\t",
                        chromCol: Int = 1,
                        skip: Int = 0
                      ): DataFrame = {

    // First, read a single line to get column names
    val sampleDf = spark.read
      .option("header", header)
      .option("delimiter", sep)
      .option("inferSchema", "true")
      .csv(file)
      .limit(1)

    // Get column names and create schema with chromosome column as String
    val columnNames = sampleDf.columns

    // Create schema ensuring chromCol is treated as String
    val schemaFields = new ListBuffer[StructField]()

    for (i <- columnNames.indices) {
      val dataType = if ((i + 1) == chromCol) StringType else sampleDf.schema(i).dataType
      schemaFields += StructField(columnNames(i), dataType, true)
    }

    val schema = StructType(schemaFields.toList)

    // Read the complete file with the custom schema
    var df = spark.read
      .format("csv")
      .option("header", header)
      .option("delimiter", sep)
      .option("mode", "DROPMALFORMED")
      .schema(schema)
      .load(file)

    // Skip rows if needed
    if (skip > 0) {
      df = df.offset(skip)
    }

    // Handle row names if needed
    if (rowNames) {
      val firstColName = df.columns(0)
      val rowNamesCol = df.columns(0)

      // Convert first column to row names (in Spark, we keep the column but can set it as an index for certain operations)
      df = df.withColumn("rowName", col(firstColName))
        .drop(firstColName)
    }

    // Replace spaces with dots in column names
    val newColumnNames = df.columns.map(colName => colName.replace(" ", "."))

    // Apply the new column names
    for (i <- df.columns.indices) {
      df = df.withColumnRenamed(df.columns(i), newColumnNames(i))
    }

    df
  }

  import java.io.File
  import scala.collection.mutable.ListBuffer

  /**
   * Function to concatenate allele counter output
   *
   * @param spark      The SparkSession instance
   * @param inputStart The prefix of the file path
   * @param inputEnd   The suffix of the file path
   * @return A DataFrame with concatenated allele count data
   */
  def concatenateAlleleCountFiles(
                                   spark: SparkSession,
                                   inputStart: String,

                                 ): DataFrame = {

    var validFiles = new ListBuffer[String]()


    // Collect valid files that exist and have data
    for (chrom <- chromosomeNames.toList) {
      val filename = s"${inputStart}${chrom}.txt"
      val file = new File(filename)

      if (file.exists() && file.length() > 0) {
        validFiles += filename
      }
    }

    if (validFiles.isEmpty) {
      return spark.emptyDataFrame
    }

    // Read and concatenate all valid files
    var resultDf: DataFrame = null

    for (file <- validFiles) {
      val currentDf = readTableGeneric(spark, file)

      if (resultDf == null) {
        resultDf = currentDf
      } else {
        resultDf = resultDf.union(currentDf)
      }
    }

    resultDf
  }

  import org.apache.spark.sql.functions.*
  import org.apache.spark.sql.{DataFrame, SparkSession}

  import java.io.File

  /**
   * Function to concatenate G1000 SNP files
   *
   * @param spark      The SparkSession instance
   * @param inputStart The prefix of the file path
   * @return A DataFrame with concatenated G1000 SNP data
   */
  def concatenateG1000SnpFiles(
                                spark: SparkSession,
                                inputStart: String,
                              ): DataFrame = {

    var resultDf: DataFrame = null


    for (chrom <- chromosomeNames.toList) {
      val filename = s"${inputStart}${chrom}.txt"
      val file = new File(filename)

      if (file.exists() && file.length() > 0) {
        // Read the file
        val chromDf = readTableGeneric(spark, filename)

        // Add chromosome column to the data
        val withChromDf = chromDf.withColumn("CHR", lit(chrom))

        // Append to the result
        if (resultDf == null) {
          resultDf = withChromDf
        } else {
          resultDf = resultDf.union(withChromDf)
        }
      }
    }
    val uniqueValues = resultDf.select("CHR").distinct()
    // Return empty DataFrame if no files were processed
    if (resultDf == null) {
      return spark.emptyDataFrame
    }

    resultDf
  }

  def createPatientDirectory(patient_name: String): Unit = {

    working_directory = s"$working_directory/$patient_name"
    plots_directory = s"$working_directory/$plots_directory"
    impute_directory = s"$working_directory/$impute_directory"
    allele_directory = s"$working_directory/$allele_directory"

    val new_directories: Seq[String] = Seq(working_directory, plots_directory, impute_directory, allele_directory)

    new_directories.foreach(checkCorrectExecution)

  }


  def setChromosomesNames(bamFile: String): Unit = {
    val command = Seq(
      "samtools",
      "view", "-H", bamFile
    ).mkString(" ")
    println(command)


    val chromosomePrefix = command.!!

    if (chromosomePrefix.contains("chr")) {
      chromosomeNames = ((1 to 22).map(n => s"chr$n") :+ "chrX").toList.par

    } else {
      chromosomeNames = ((1 to 22) :+ "X").toList.par
    }

  }

  def isMale(sampleName: String): Boolean = {
    val spark = SparkSession.builder()
      .appName("Battenberg")
      .master("local[*]")
      .getOrCreate()

    try {

      val samtoolsCmd = Seq("samtools", "idxstats", sampleName)
      val sortCmd = Seq("sort")
      val fullCommand = (samtoolsCmd #| sortCmd)
      val output = Try(fullCommand.!!).getOrElse {
        spark.stop()
        throw new RuntimeException(s"Failed to execute command: $fullCommand")
      }


      val schema = StructType(Array(
        StructField("Chromosome", StringType, nullable = false),
        StructField("Length", LongType, nullable = false),
        StructField("Mapped", LongType, nullable = false),
        StructField("Unmapped", LongType, nullable = false)
      ))

      // Convert output to Row objects
      val rows = output.trim.split("\n").map { line =>
        val fields = line.split("\t")
        if (fields.length == 4) {
          Row(fields(0), fields(1).toLong, fields(2).toLong, fields(3).toLong)
        } else {
          Row("invalid", 0L, 0L, 0L) // Handle invalid lines
        }
      }


      // Create DataFrame from rows and schema
      val df = spark.createDataFrame(
        spark.sparkContext.parallelize(rows.toSeq),
        schema
      )

      val filteredDf = df.filter(!col("Chromosome").contains("_") && col("Chromosome") =!= "chrM")
        .withColumn("Length_per_Read", col("Mapped") / col("Length"))


      val chrXStatsRows = filteredDf.filter(col("Chromosome") === "chrX")
        .select("Length_per_Read")
        .collect()

      if (chrXStatsRows.isEmpty) {
        spark.stop()
        throw new RuntimeException("Chromosome X not found in the data")
      }

      val chrXLengthPerRead = chrXStatsRows(0).getDouble(0)

      val statsRow = filteredDf.agg(
        sum("Length_per_Read").as("totalLengthPerRead"),
        count("*").as("count")
      ).collect()(0)

      val totalLengthPerRead = statsRow.getDouble(0)
      val rowCount = statsRow.getLong(1)

      val avgLengthPerRead = (totalLengthPerRead - chrXLengthPerRead) / (rowCount - 1)

      // Compare values
      val comparison = Math.abs(chrXLengthPerRead - avgLengthPerRead)
      val result = comparison < 0.1

      result

    } catch {
      case e: Exception =>
        spark.stop()
        throw e
    }
  }


  private def checkCorrectExecution(path_to_directory: String): Unit = {
    val dir = new File(path_to_directory)
    if (!dir.exists()) {
      val created = dir.mkdir()
      require(created, s"Patient directory couldn't be created")
    }
  }








