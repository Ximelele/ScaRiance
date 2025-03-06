import scala.sys.process.*
//import scala.collection.parallel.CollectionConverters._
case class PrepareWgs():

//  def prepareWgs(utils: Utils): Unit={
//
//
//  }
  def getAlleleCounts(bamFile: String, outputFile: String, g1000Loci: String, minBaseQual: Int = 20, minMapQual: Int = 35, alleleCounter: String = "alleleCounter"): Unit = {
    var alleleCounterCommand = Seq(
      alleleCounter,
      "-b", bamFile,
      "-l", g1000Loci,
      "-o", outputFile,
      "-m", minBaseQual.toString,
      "-q", minMapQual.toString
    ).mkString(" ")

//    val counterVersion: String = s"$alleleCounter --version".!!.trim
//    println(s"AlleleCounter version: $counterVersion")
//
//    if counterVersion.substring(0, 1).toInt >= 4 then
//      alleleCounterCommand = s"$alleleCounterCommand --dense-snps"
//
//    println(s"Executing command: $alleleCounterCommand")
//
//    val exitCode: Int = alleleCounterCommand.!
//
//    require(exitCode == 0, s"Command failed with exit code $exitCode")
  }

  def test(): Unit = {
    println("Hehe")
  }