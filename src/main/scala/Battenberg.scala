import org.apache.spark.sql.SparkSession

import scala.collection.parallel.CollectionConverters.*

case class Battenberg(control_file: String, tumour_file: String):

  private val utils = Utils()
  private var male: Boolean = true
  private val prepare_Wgs = PrepareWgs()

  private def setDefaultValues(): Unit = {
    this.utils.setChromosomesNames(control_file)
    this.utils.createPatientDirectory(control_file.split("/").last)
  }

  def run(): Unit = {
    this.setDefaultValues()
    val spark = SparkSession.builder()
      .appName("Battenberg")
      .master("local[*]")
      .getOrCreate()

    this.male = this.utils.isMale(tumour_file)
    println("Idem riesit allele")
    prepare_Wgs.prepareWgs(utils = utils, controlFile = control_file, tumourFile = tumour_file)

    spark.stop()

  }


