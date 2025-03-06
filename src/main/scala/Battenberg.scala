import scala.collection.parallel.CollectionConverters.*

case class Battenberg(control_file: String, tumour_file: String):

  private val utils = Utils()


  def isMale: Boolean = {

    true
  }

  private def haf(chromosome: Any): Unit = {
    println(s"Processing $chromosome on thread ${Thread.currentThread().getName}")
  }

  private def setDefaultValues(): Unit = {
    this.utils.setChromosomesNames(control_file)
    this.utils.createPatientDirectory(control_file.split("/").last)

    this.utils.chromosomeNames.foreach(chromsome => {

      haf(chromsome)
    }
    )

  }

  def run(): Unit = {
    this.setDefaultValues()

  }


