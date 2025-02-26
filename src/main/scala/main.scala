
@main
def main(): Unit =

  val mixedSeq: Seq[Any] = (1 to 22).map(n => s"chr$n") :+ "chrX"


  for (chromosome <- mixedSeq) {
    println(chromosome)
  }

  val battenberg = Battenberg()

  battenberg.getChromosomesNames("/data/I002.007.WGS.C.bam")
//  val prapareWgs = PrepareWgs()
//  val loci = "/app/references38/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chr"
//
//
//  for (chromosome_name <- 1 to 10)
//
//    prapareWgs.getAlleleCounts(bamFile = "/data/I002.007.WGS.C.bam", outputFile = s"/app/_alleleFrequencies_chr$chromosome_name.txt", g1000Loci = s"${loci}$chromosome_name.txt")
