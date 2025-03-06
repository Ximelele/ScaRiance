import java.io.File
import scala.collection.parallel.CollectionConverters.*
import scala.collection.parallel.ParSeq
import scala.sys.process.*

case class Utils():
  var working_directory: String = sys.props("user.dir")
  var impute_directory: String = "Impute"
  var plots_directory: String = "Plots"
  var allele_directory: String = "Allele_frequencies"
  var chromosomeNames: ParSeq[Any] = ParSeq()

  def isMale: Boolean = {
    true
  }


  private def checkCorrectExecution(path_to_directory: String): Unit = {
    val dir = new File(path_to_directory)
    if (!dir.exists()) {
      val created = dir.mkdir()
      require(created, s"Patient directory couldn't be created")
    }
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
    //    val command = Seq(
    //      "samtools",
    //      "view", "-H", bamFile
    //    ).mkString(" ")
    //    println(command)
    //
    //
    //    val chromosomePrefix = command.!!
    //    val chromosomePrefix = "chr"

    //    if (chromosomePrefix.contains("chr")) {
    chromosomeNames = ((1 to 22).map(n => s"chr$n") :+ "chrX").toList.par

    //    } else {
    //      chromosomeNames = (1 to 22) :+ "chrX"
    //    }
  }






