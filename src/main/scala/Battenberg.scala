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
    this.utils.chromosomeNames.foreach(chrom => {
      println(chrom)
      //impute.runHaplotyping(spark, chrom, this.utils)
      haplotype.getChromosomeBafs(spark = spark, SNP_file = s"${utils.allele_directory}/${utils.tumourName}_alleleFrequencies_$chrom.txt", haplotype_File = s"${utils.impute_directory}/${utils.tumourName}_impute_output_${chrom}_allHaplotypeInfo.txt", utils = utils, output_file = s"${utils.impute_directory}/${utils.tumourName}_impute_output_${chrom}_heterozygousMutBAFs_haplotyped.txt", minCounts = 10)
    })


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


