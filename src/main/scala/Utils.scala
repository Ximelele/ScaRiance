import org.apache.spark.sql.*
import org.apache.spark.sql.functions.*
import org.apache.spark.sql.types.*

import java.io.File
import java.nio.file.{Files, Paths, StandardCopyOption}
import scala.collection.mutable.ListBuffer
import scala.collection.parallel.CollectionConverters.*
import scala.collection.parallel.ParSeq
import scala.sys.process.*
import scala.util.{Random, Try}

case class Utils():
  val referenciesFile: References = References()
  var working_directory: String = sys.props("user.dir")
  var impute_directory: String = "Impute"
  var plots_directory: String = "Plots"
  var allele_directory: String = "Allele_frequencies"
  var chromosomeNames: ParSeq[Any] = ParSeq()
  var tumourName: String = ""
  var controlName: String = ""
  var is_male: Boolean = true


  def saveSingleFile(df: DataFrame, outputPath: String, header: Boolean = true): Unit = {

    val tempDir = s"${outputPath}_temp"
    df.coalesce(1).write
      .option("header", header)
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


  /**
   * Generic reading function tailored for reading in genomic data using Spark
   *
   * @param spark    The SparkSession instance
   * @param file     Filename of the file to read in
   * @param header   Whether the file contains a header (Default: true)
   * @param rowNames Whether the file contains row names (Default: false)
   * @param sep      Column separator (Default: "\t")
   * @return A DataFrame with contents of the file
   */
  private def readTableGeneric(
                                spark: SparkSession,
                                file: String,
                                schema: StructType
                              ): DataFrame = {

    var df = spark.read
      .format("csv")
      .option("header", true)
      .option("delimiter", "\t")
      .option("mode", "DROPMALFORMED")
      .schema(schema)
      .load(file)


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
   * @return A DataFrame with concatenated allele count data
   */
  def concatenateAlleleCountFiles(
                                   spark: SparkSession,
                                   inputStart: String,

                                 ): DataFrame = {

    val validFiles = new ListBuffer[String]()

    val alleleCountSchema = StructType(Seq(
      StructField("#CHR", StringType, nullable = true),
      StructField("POS", IntegerType, nullable = true),
      StructField("Count_A", IntegerType, nullable = true),
      StructField("Count_C", IntegerType, nullable = true),
      StructField("Count_G", IntegerType, nullable = true),
      StructField("Count_T", IntegerType, nullable = true),
      StructField("Good_depth", IntegerType, nullable = true)
    ))

    // Collect valid files that exist and have data
    for (chrom <- chromosomeNames.toList) {
      val filename = s"$inputStart$chrom.txt"
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
      val currentDf = readTableGeneric(spark, file, alleleCountSchema)

      if (resultDf == null) {
        resultDf = currentDf
      } else {
        resultDf = resultDf.union(currentDf)
      }
    }
    resultDf.show()
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


    val g1000SnpSchema = StructType(Seq(
      StructField("position", IntegerType, nullable = true),
      StructField("a0", IntegerType, nullable = true),
      StructField("a1", IntegerType, nullable = true)
    ))

    for (chrom <- chromosomeNames.toList) {
      val filename = s"$inputStart$chrom.txt"
      val file = new File(filename)

      if (file.exists() && file.length() > 0) {
        // Read the file
        val chromDf = readTableGeneric(spark, filename, g1000SnpSchema)

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
    resultDf.show()
    resultDf
  }


  def concatenateBAFfiles(spark: SparkSession, inputStart: String): Unit = {

    var resultDf: DataFrame = null

    val bafSegment = StructType(Seq(
      StructField("Chromosome", StringType, nullable = true),
      StructField("Position", IntegerType, nullable = true),
      StructField(tumourName, DoubleType, nullable = true)
    ))


    for (chrom <- chromosomeNames.toList) {
      val filename = s"$inputStart${chrom}_heterozygousMutBAFs_haplotyped.txt"

      val file = new File(filename)

      if (file.exists() && file.length() > 0) {

        val chromDf = readTableGeneric(spark, filename,bafSegment)

        // Append to the result
        if (resultDf == null) {
          resultDf = chromDf
        } else {
          resultDf = resultDf.union(chromDf)
        }
      }
    }
    val extractChr = regexp_extract(col("Chromosome"), "^chr?(\\d+|X)", 1)

    resultDf = resultDf.orderBy(
      when(extractChr === "X", 100).otherwise(extractChr.cast("int")),
      col("Position").cast("int")
    )
    resultDf.show()
    saveSingleFile(resultDf, s"$impute_directory/${tumourName}_heterozygousMutBAFs_haplotyped.txt")
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
    //
    if (chromosomePrefix.contains("chr")) {
      chromosomeNames = ((1 to 22).map(n => s"chr$n") :+ "chrX").toList.par
      //
    } else {
      chromosomeNames = ((1 to 22) :+ "X").toList.par
    }

    println(chromosomeNames)

  }


  def isMale(sampleName: String, spark: SparkSession): Boolean = {

    try {

      val samtoolsCmd = Seq("samtools", "idxstats", sampleName)
      val sortCmd = Seq("sort")
      val fullCommand = samtoolsCmd #| sortCmd
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

      println("Patient is male: ",result)

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








