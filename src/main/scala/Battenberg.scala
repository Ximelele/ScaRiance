import org.apache.spark.sql.SparkSession

import scala.collection.parallel.CollectionConverters.*

case class Battenberg(control_file: String, tumour_file: String):

  private val utils = Utils()
  private val prepare_Wgs = PrepareWgs()
  private var male: Boolean = true

  def run(): Unit = {
    this.setDefaultValues()
    val spark = SparkSession.builder()
      .appName("Battenberg")
      .master("local[*]")
      .getOrCreate()

    this.male = this.utils.isMale(tumour_file)
    concatenateAlleleCountFiles(inputStart = s"${this.utils.allele_directory}/${controlFile.split("/").last}_alleleFrequencies_", inputEnd = ".txt", chrNames = this.utils.chromosomeNames)
    //    prepare_Wgs.prepareWgs(utils = utils, controlFile = control_file, tumourFile = tumour_file)

    spark.stop()

  }

  private def setDefaultValues(): Unit = {
    this.utils.setChromosomesNames(control_file)
    this.utils.createPatientDirectory(control_file.split("/").last)
  }


