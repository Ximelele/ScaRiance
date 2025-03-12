import org.apache.spark.sql.functions.*
import org.apache.spark.sql.types.*
import org.apache.spark.sql.{DataFrame, Row, SparkSession}

import java.io.File
import scala.collection.mutable.ArrayBuffer
import scala.collection.parallel.CollectionConverters.*
import scala.collection.parallel.ParSeq
import scala.io.Source
import scala.sys.process.*
import scala.util.Try

case class Utils():
  val g1000prefix: String = "/app/references38/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_"
  var working_directory: String = sys.props("user.dir")
  var impute_directory: String = "Impute"
  var plots_directory: String = "Plots"
  var allele_directory: String = "Allele_frequencies"
  var chromosomeNames: ParSeq[Any] = ParSeq()


  def concatenateAlleleCountFiles(inputStart: String, inputEnd: String, chrNames: List[String]): DataFrame = {
    import org.apache.spark.sql.{DataFrame, SparkSession}

    val spark = SparkSession.builder()
      .appName("Battenberg")
      .getOrCreate()

    // Find valid input files
    val infiles = chrNames.flatMap { chrom =>
      val filename = s"${inputStart}${chrom}${inputEnd}"
      val file = new File(filename)
      if (file.exists() && file.length() > 0) {
        Some(filename)
      } else {
        None
      }
    }

    // Read and combine all files
    if (infiles.isEmpty) {
      // Return empty DataFrame with appropriate schema
      spark.emptyDataFrame
    } else {
      // Read first file to get schema
      val firstDf = readTableGeneric(infiles.head)

      // Read remaining files and union them
      val allDfs = infiles.tail.foldLeft(firstDf) { (combinedDf, file) =>
        val df = readTableGeneric(file)
        combinedDf.union(df)
      }

      allDfs.allDfs.write
        .mode("overwrite") // Can also be "append"
        .option("delimiter", "\t") // Set delimiter to tab character
        .option("header", "true") // Optionally include header row
        .csv("test.txt")
    }
  }

  def readTableGeneric(
                        file: String,
                        header: Boolean = true,
                        rowNames: Boolean = false,
                        stringsAsFactor: Boolean = false, // kept for legacy purposes but not used
                        sep: String = "\t",
                        chromCol: List[Int] = List(1),
                        skip: Int = 0
                      ): DataFrame = {

    // Skip initial lines if needed
    val iterator = Source.fromFile(file).getLines().drop(skip)

    // Read first line to get header
    val firstLine = if (iterator.hasNext) iterator.next() else ""
    val headerColumns = if (header) firstLine.split(sep).map(_.replaceAll(" ", ".")) else
      (1 to firstLine.split(sep).length).map(i => s"V$i").toArray

    // Read the rest of the file
    val rows = ArrayBuffer[Array[String]]()

    // Add first line to data if it's not a header
    if (!header) {
      rows += firstLine.split(sep)
    }

    // Read remaining lines
    while (iterator.hasNext) {
      val line = iterator.next()
      rows += line.split(sep)
    }


    // Create SparkSession
    val spark = SparkSession.builder()
      .appName("Battenberg")
      .getOrCreate()

    // Infer schema (with special handling for chromosome columns)
    val schemaFields = headerColumns.zipWithIndex.map { case (colName, idx) =>
      if (chromCol.contains(idx + 1)) {
        StructField(colName, StringType, nullable = true)
      } else {
        // Simple type inference - could be enhanced
        val sampleValues = rows.take(10).map(row => if (row.length > idx) row(idx) else "")
        val isNumeric = sampleValues.forall(s => s.isEmpty || s.matches("^-?\\d+(\\.\\d+)?$"))
        if (isNumeric) StructField(colName, DoubleType, nullable = true)
        else StructField(colName, StringType, nullable = true)
      }
    }

    val schema = StructType(schemaFields)

    // Create rows
    val rowsRDD = spark.sparkContext.parallelize(rows.map(arr => Row(arr: _*)))

    // Create DataFrame
    var df = spark.createDataFrame(rowsRDD, schema)

    // Handle row names if needed
    if (rowNames) {
      // Create a new column with row names
      import org.apache.spark.sql.functions.*
      val rowNamesColumn = df.columns(0)
      df = df.drop(rowNamesColumn)
      // Note: In Spark DataFrames, there's no direct equivalent to R's row names
    }

    df
  }

  def createPatientDirectory(patient_name: String): Unit = {

    working_directory = s"$working_directory/$patient_name"
    plots_directory = s"$working_directory/$plots_directory"
    impute_directory = s"$working_directory/$impute_directory"
    allele_directory = s"$working_directory/$allele_directory"

    val new_directories: Seq[String] = Seq(working_directory, plots_directory, impute_directory, allele_directory)

    new_directories.foreach(checkCorrectExecution)

  }

  // Helper function to create a list of chromosomes (1-22 + X)
  def createChromosomeList(): Seq[String] = {
    val autosomes = (1 to 22).map(i => s"chr$i")
    autosomes :+ "chrX"
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








