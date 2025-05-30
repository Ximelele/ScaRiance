import org.apache.spark.sql.functions.col
import org.apache.spark.sql.{DataFrame, SparkSession}
import org.apache.spark.storage.StorageLevel


case class Impute():

  private def parseImputeFile(spark: SparkSession, is_male: Boolean, chrom: String, impute_file: String): DataFrame = {
    var impute_info = spark.read.option("header", "false")
      .option("delimiter", "\t")
      .csv(impute_file)

    impute_info = impute_info.toDF("chrom", "impute_legend", "genetic_map", "impute_hap", "start", "end", "is_par")

    if (is_male) {
      impute_info = impute_info.filter(col("is_par") === 1)
    }

    impute_info = impute_info.filter(col("chrom") === chrom)
    impute_info
  }


  private def writeBeagleAsImpute(spark: SparkSession, beagle_file: String, output_file: String, utils: Utils): Unit = {
    var beagle_output = spark.read.option("header", "false").option("comment", "#").option("delimiter", "\t").csv(beagle_file)
    import org.apache.spark.sql.functions.*
    import org.apache.spark.sql.types.*

    val customSchema = StructType(Array(
      StructField("CHROM", StringType, true),
      StructField("POS", IntegerType, true),
      StructField("ID", StringType, true),
      StructField("REF", StringType, true),
      StructField("ALT", StringType, true),
      StructField("QUAL", DoubleType, true),
      StructField("FILTER", StringType, true),
      StructField("INFO", StringType, true),
      StructField("FORMAT", StringType, true),
      StructField("SAMP0001", StringType, true)
    ))
    beagle_output = beagle_output.toDF(customSchema.fieldNames: _*)
      .select(
        col("CHROM"),
        col("POS").cast(IntegerType),
        col("ID"),
        col("REF"),
        col("ALT"),
        col("QUAL").cast(DoubleType),
        col("FILTER"),
        col("INFO"),
        col("FORMAT"),
        col("SAMP0001")
      )

    beagle_output = beagle_output.withColumn("genotype_split", split(col("SAMP0001"), "|"))
      .withColumn("allele1", col("genotype_split").getItem(0).cast("integer"))
      .withColumn("allele2", col("genotype_split").getItem(2).cast("integer"))
      .drop("genotype_split")

    beagle_output = beagle_output.withColumn("row_num", monotonically_increasing_id() + 1)

    beagle_output = beagle_output
      .withColumn("snp_index", concat(lit("snp_index"), col("row_num").cast("string")))
      .withColumn("rs_index", concat(lit("rs_index"), col("row_num").cast("string")))

    beagle_output = beagle_output.select(col("snp_index"), col("rs_index"), col("POS"), col("REF"), col("ALT"), col("allele1"), col("allele2"))

    utils.saveSingleFile(beagle_output, output_file, false)
  }

  private def generateImputeInput(spark: SparkSession, chromosome: String, utils: Utils): String = {

    val tumourAlleleCountsFile = s"${utils.allele_directory}/${utils.tumourName}_alleleFrequencies_$chromosome.txt"
    val normalAlleleCountsFile = s"${utils.allele_directory}/${utils.controlName}_alleleFrequencies_$chromosome.txt"
    val impute_info = parseImputeFile(spark, utils.is_male, chromosome, utils.referenciesFile.impute_file)
    import spark.implicits.*
    val impute_legend_files = impute_info.select("impute_legend").as[String].collect()

    var known_SNPs = spark.read.option("header", "true").option("sep", " ").csv(impute_legend_files(0))

    if (impute_legend_files.length > 1) {
      for (idx <- 1 until impute_legend_files.length) {
        val next_SNP = spark.read.option("header", "true").option("sep", " ").csv(impute_legend_files(idx))
        known_SNPs = known_SNPs.union(next_SNP)
      }
    }

    var problemSNPs = spark.read.option("header", "true").option("sep", "\t").csv(utils.referenciesFile.problemLociFile)


    problemSNPs = problemSNPs.filter(col("Chr") === chromosome).select("Pos")
    known_SNPs = known_SNPs.join(problemSNPs, known_SNPs("position") === problemSNPs("Pos"), "left_anti")

    problemSNPs.unpersist()

    var snp_data = spark.read.option("header", "true").option("sep", "\t").csv(tumourAlleleCountsFile)
    var normal_snp_data = spark.read.option("header", "true").option("sep", "\t").csv(normalAlleleCountsFile)

    snp_data = snp_data.columns.filter(_ != "POS").foldLeft(snp_data) {
      case (df, colName) => df.withColumnRenamed(colName, s"snp_$colName")
    }

    normal_snp_data = normal_snp_data.columns.filter(_ != "POS").foldLeft(normal_snp_data) {
      case (df, colName) => df.withColumnRenamed(colName, s"normal_$colName")
    }

    snp_data = snp_data.join(normal_snp_data, Seq("POS"), "inner")

    snp_data = snp_data.join(
      known_SNPs.select("position", "id"),
      snp_data("POS") === known_SNPs("position"),
      "inner"
    ).orderBy(col("POS").cast("int"))

    import org.apache.spark.sql.functions.*


    var withBAF = snp_data
      .join(known_SNPs.select("position", "a0", "a1"),
        snp_data("position") === known_SNPs("position"),
        "inner")
      .withColumn("ref_count",
        when(col("a0") === "A", col("snp_Count_A"))
          .when(col("a0") === "C", col("snp_Count_C"))
          .when(col("a0") === "G", col("snp_Count_G"))
          .when(col("a0") === "T", col("snp_Count_T"))
          .otherwise(lit(0)).cast("double"))
      .withColumn("alt_count",
        when(col("a1") === "A", col("snp_Count_A"))
          .when(col("a1") === "C", col("snp_Count_C"))
          .when(col("a1") === "G", col("snp_Count_G"))
          .when(col("a1") === "T", col("snp_Count_T"))
          .otherwise(lit(0)).cast("double"))
      .withColumn("BAF",
        when(col("alt_count") + col("ref_count") > 0,
          col("alt_count") / (col("alt_count") + col("ref_count")))
          .otherwise(lit(0.0)))


    withBAF = withBAF.drop("position")

    // todo create class or enum to set this
    val heterozygousFilter = 0.1
    val minBaf = math.min(heterozygousFilter, 1.0 - heterozygousFilter)
    val maxBaf = math.max(heterozygousFilter, 1.0 - heterozygousFilter)


    withBAF = withBAF
      .withColumn("genotype_0",
        when(col("BAF") <= minBaf, 1).otherwise(0))
      .withColumn("genotype_1",
        when(col("BAF") > minBaf && col("BAF") < maxBaf, 1).otherwise(0))
      .withColumn("genotype_2",
        when(col("BAF") >= maxBaf, 1).otherwise(0))

    withBAF = withBAF.withColumn("ID", lit("."))
      .withColumn("QUAL", lit("."))
      .withColumn("FILTER", lit("PASS"))
      .withColumn("INFO", lit("."))
      .withColumn("FORMAT", lit("GT"))
      .withColumn("SAMP0001", concat_ws("-", col("genotype_0"), col("genotype_1"), col("genotype_2")))
      .withColumn("SAMP0001",
        when(col("SAMP0001") === "1-0-0", "0/0")
          .when(col("SAMP0001") === "0-1-0", "0/1")
          .when(col("SAMP0001") === "0-0-1", "1/1")
      )

    withBAF = withBAF.filter(!(col("SAMP0001") === "0-0-0"))

    withBAF = withBAF
      .select(
        col("snp_#CHR").as("CHROM"),
        col("POS"),
        col("ID"),
        col("a0").as("REF"),
        col("a1").as("ALT"),
        col("QUAL"),
        col("FILTER"),
        col("INFO"),
        col("FORMAT"),
        col("SAMP0001"),
      ).orderBy(col("POS").cast("int"))


    val outputFile = s"${utils.impute_directory}/${utils.tumourName}_beagle5_input_$chromosome.txt"
    withBAF = withBAF.persist(StorageLevel.MEMORY_AND_DISK)
    utils.saveSingleFile(withBAF, outputPath = outputFile)

    val headerContent = Seq(
      "##fileformat=VCFv4.2",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "##reference=hg38"
    )
    import scala.sys.process.*

    val cmd = Seq("bash", "-c",
      s"(echo '${headerContent.mkString("\n")}'; cat $outputFile) > $outputFile.temp && mv $outputFile.temp $outputFile")

    cmd.!

    val chrom_cmd = s"sed -i 's/^CHROM/#CHROM/g' $outputFile"
    chrom_cmd.!

    val chrom_prefix = s"sed -i 's/^chr//g' $outputFile"

    chrom_prefix.!

    outputFile
  }

  private def runBeagle5(vcf_path: String, utils: Utils, output_path: String, beaglemaxmem: Int = 10,
                         beaglewindow: Int = 40,
                         beagleoverlap: Int = 4, nthreads: Int = 1, window: Int = 1, overlap: Int = 1, chromosome: String): Unit = {


    val cmd = Seq("java",
      s"-Xmx${beaglemaxmem}g",
      s"-Xms${beaglemaxmem}g",
      s"-XX:+UseParallelGC",
      s"-jar ${utils.referenciesFile.beaglejar}",
      s"gt=$vcf_path",
      s"ref=${utils.referenciesFile.beagleref.replace("CHROMNAME", chromosome)}",
      s"out=$output_path",
      s"map=${utils.referenciesFile.beagleplink.replace("CHROMNAME", chromosome)}",
      s"nthreads=$nthreads",
      s"window=$beaglewindow",
      s"overlap=$beagleoverlap",
      "impute=false").mkString(" ")


    import scala.sys.process.*

    cmd.!
  }

  def runHaplotyping(spark: SparkSession, chromosome: String, utils: Utils): Unit = {


    val vcfbeagle_path = generateImputeInput(spark, chromosome, utils)
    val outbeagle_path = s"${utils.impute_directory}/${utils.tumourName}_beagle5_output_chr${chromosome.stripPrefix("chr")}.txt"


    runBeagle5(vcf_path = vcfbeagle_path, utils = utils, output_path = outbeagle_path, chromosome = chromosome)


    val outfile = s"${utils.impute_directory}/${utils.tumourName}_impute_output_chr${chromosome.stripPrefix("chr")}_allHaplotypeInfo.txt"

    val vcfout = s"$outbeagle_path.vcf.gz"

    writeBeagleAsImpute(spark, vcfout, outfile, utils)


  }




