import scala.sys.process.*

case class Battenberg():

  //  val chromosomeNames : Seq[Int|String]
  def isMale: Boolean = {

    return true
  }

  def getChromosomesNames(bamFile: String): Unit = {
    val command = s"samtools view -H $bamFile | grep '^@SQ' | cut -f2 | sed 's/SN://' | head -1"
    println(command)

    val chromosomePrefix = command.!!
    println(chromosomePrefix)
  }
