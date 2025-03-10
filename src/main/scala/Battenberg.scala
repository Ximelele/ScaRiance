import org.apache.spark.sql.SparkSession

import scala.collection.parallel.CollectionConverters.*

case class Battenberg(control_file: String, tumour_file: String):

  private val utils = Utils()


  private def setDefaultValues(): Unit = {
    this.utils.setChromosomesNames(control_file)
    this.utils.createPatientDirectory(control_file.split("/").last)

    //    this.utils.chromosomeNames.foreach(chromsome => {
    //
    //      haf(chromsome)
    //    }
    //    )

  }

  def run(): Unit = {
    this.setDefaultValues()
    val spark = SparkSession.builder()
      .appName("Battenberg")
      .master("local[*]")
      .getOrCreate()

    this.utils.isMale(tumour_file)

  }


