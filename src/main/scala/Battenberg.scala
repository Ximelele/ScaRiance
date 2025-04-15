import org.apache.spark.sql.SparkSession

import scala.collection.parallel.CollectionConverters.*

case class Battenberg(control_file: String, tumour_file: String):

  private val utils = Utils()
  private val prepare_Wgs = PrepareWgs()
  private var male: Boolean = true

  def run(): Unit = {
    this.setDefaultValues()


    //    this.male = this.utils.isMale(tumour_file)
    prepare_Wgs.prepareWgs(utils = utils, controlFile = control_file, tumourFile = tumour_file)
    println("Macka Macka")
    spark.stop()

  }

  private def setDefaultValues(): Unit = {
    this.utils.tumourName = tumour_file.split("/").last
    this.utils.controlName = control_file.split("/").last
    this.utils.setChromosomesNames(control_file)
    this.utils.createPatientDirectory(this.utils.controlName)
  }


