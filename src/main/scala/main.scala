
@main
def main(): Unit =
  val prapareWgs = Prepare_wgs()
  val loci = "/app/references38/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chr"
  for (chromosome_name <- 1 to 2)

    prapareWgs.getAlleleCounts(bamFile = "/data/I002.007.WGS.C.bam", outputFile = "./", g1000Loci = s"${loci}$chromosome_name")
  println("Mackaaaa")

