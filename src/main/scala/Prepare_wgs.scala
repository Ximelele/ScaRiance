import scala.sys.process.*

class Prepare_wgs:
  def getAlleleCounts(bamFile: String, outputFile: String, g1000Loci: String, minBaseQual: Int, minMapQual: Int, alleleCounter: String): Unit = {
    var alleleCounterCommand = Seq(
      alleleCounter,
      "-b", bamFile,
      "-l", g1000Loci,
      "-o", outputFile,
      "-m", minBaseQual.toString,
      "-q", minMapQual.toString
    ).mkString(" ")

    val counterVersion: String = s"$alleleCounter --version".!!.trim
    println(s"AlleleCounter version: $counterVersion")

    if counterVersion.substring(0, 1).toInt >= 4 then
      alleleCounterCommand = s"$alleleCounterCommand --dense-snps"

    println(s"Executing command: $alleleCounterCommand")

    val exitCode: Int = alleleCounterCommand.!
    // Stop the program if the exit code is not 0.
    require(exitCode == 0, s"Command failed with exit code $exitCode")
  }