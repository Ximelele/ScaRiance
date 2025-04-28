
import org.apache.spark.sql.functions.{col, lit, when}
import org.apache.spark.sql.{Column, DataFrame, SparkSession}

case class Haplotype():


  def getChromosomeBafs(spark: SparkSession, SNP_file: String, haplotype_File: String, utils: Utils, output_file: String, minCounts: Int = 1): Unit = {

    var snp_data = spark.read.option("header", "true").option("sep", "\t").csv(SNP_file)
    var variant_data = spark.read.option("header", "false").option("sep", "\t").csv(haplotype_File)

    variant_data = variant_data.toDF("id1", "id2", "POS_var", "REF", "ALT", "allele1", "allele2")

    variant_data = variant_data.filter(col("allele1") =!= col("allele2"))

    variant_data.show()


    snp_data = snp_data.join(variant_data, snp_data("POS") === variant_data("POS_var"), "inner")
    import org.apache.spark.sql.types._
    val schema = StructType(Array(
      StructField("Chromosome", StringType, nullable = true),
      StructField("Position", IntegerType, nullable = true),
      StructField(utils.tumourName, StringType, nullable = true) // Using variable for column name
    ))
    val emptyDF: DataFrame = spark.createDataFrame(
      spark.sparkContext.emptyRDD[org.apache.spark.sql.Row],
      schema
    )
    if (snp_data.count() == 0) {
      utils.saveSingleFile(emptyDF, output_file)
      return
    }
    snp_data = snp_data.withColumn("ref_count",
        when(col("REF") === "A", col("Count_A"))
          .when(col("REF") === "C", col("Count_C"))
          .when(col("REF") === "G", col("Count_G"))
          .when(col("REF") === "T", col("Count_T"))
          .otherwise(lit(0)).cast("double"))
      .withColumn("alt_count",
        when(col("ALT") === "A", col("Count_A"))
          .when(col("ALT") === "C", col("Count_C"))
          .when(col("ALT") === "G", col("Count_G"))
          .when(col("ALT") === "T", col("Count_T"))
          .otherwise(lit(0)).cast("double"))
      .withColumn("denom", col("ref_count") + col("alt_count"))
      .filter(col("denom") >= minCounts)


    if (snp_data.count() == 0) {
      utils.saveSingleFile(emptyDF, output_file)
      return
    }

    snp_data = snp_data.withColumn("hetMutBafs", col("alt_count") / col("denom"))

    snp_data.show()
    snp_data = snp_data.select(col("#CHR").as("Chromosome"), col("POS").as("Position"), col("hetMutBafs").as(utils.tumourName))

    utils.saveSingleFile(snp_data, output_file)

  }
