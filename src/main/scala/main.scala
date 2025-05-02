
object Main {
  def main(args: Array[String]): Unit = {

    var controlFile: Option[String] = None
    var tumourFile: Option[String] = None
    var skipAlleleCounting: Boolean = false
    var skipImputation: Boolean = false
    var skipSegmentation: Boolean = false

    // Parse arguments
    var i = 0
    while (i < args.length) {
      args(i) match {
        case "-control_file" =>
          if (i + 1 < args.length) {
            controlFile = Some(args(i + 1))
            i += 2
          } else {
            println("Error: -control_file requires a path value")
            printUsage()
            sys.exit(1)
          }

        case "-tumour_file" =>
          if (i + 1 < args.length) {
            tumourFile = Some(args(i + 1))
            i += 2
          } else {
            println("Error: -tumour_file requires a path value")
            printUsage()
            sys.exit(1)
          }

        case "-skip_allele_counting" =>
          skipAlleleCounting = true
          i += 1

        case "-skip_imputation" =>
          skipImputation = true
          i += 1

        case "-skip_segmentation" =>
          skipSegmentation = true
          i += 1

        case unknown =>
          println(s"Unknown option: $unknown")
          printUsage()
          sys.exit(1)
      }
    }
    import scala.reflect.io.File

    def fileExists(path: String): Boolean = {
      File(path).exists
    }
    // Check for required arguments
    if (controlFile.isEmpty) {
      println("Error: Missing required -control_file argument")
      printUsage()
      sys.exit(1)
    }
    if (!fileExists(controlFile.get)) {
      println("Error: Control file path doesnt exists")
      sys.exit(1)
    }

    if (tumourFile.isEmpty) {
      println("Error: Missing required -tumour_file argument")
      printUsage()
      sys.exit(1)
    }
    if (!fileExists(tumourFile.get)) {
      println("Error: Control file path doesnt exists")
      sys.exit(1)
    }

    // Create ScaRiance instance with all parameters
    val scaRiance = ScaRiance(
      control_file = controlFile.get,
      tumour_file = tumourFile.get,
      skip_allele_counting = skipAlleleCounting,
      skip_imputation = skipImputation,
      skip_segmentation = skipSegmentation
    )


    scaRiance.run()

    def printUsage(): Unit = {
      println("Usage: ScaRiance [options]")
      println("Required options:")
      println("  -control_file <path>      Path to control BAM file")
      println("  -tumour_file <path>       Path to tumour BAM file")
      println("Optional switches:")
      println("  -skip_allele_counting     Skip the allele counting step")
      println("  -skip_imputation          Skip the imputation step")
      println("  -skip_segmentation        Skip the segmentation step")
    }
  }
}