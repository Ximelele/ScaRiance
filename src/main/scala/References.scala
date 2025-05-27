class References:
  var g1000prefix = ""
  var g1000alleleprefix = ""
  var impute_file: String = ""
  var problemLociFile: String = ""
  var beagleref: String = ""
  var beagleplink: String = ""
  var beaglejar: String = ""


  def setPostfixDirectory(): Unit = {
    g1000prefix = "/app/references38/1000G_loci_hg38_chr/1kg.phase3.v5a_GRCh38nounref_loci_"
    g1000alleleprefix = "/app/references38/1000G_loci_hg38_chr/1kg.phase3.v5a_GRCh38nounref_allele_index_"
    impute_file = "/app/references38/imputation_chr/impute_info.txt"
    problemLociFile = "/app/references38/probloci_chr/probloci.txt.gz"
    beagleref = "/app/references38/beagle_chr/CHROMNAME.1kg.phase3.v5a_GRCh38nounref.vcf.gz"
    beaglejar = "/app/references38/beagle_chr/beagle.08Feb22.fa4.jar"
    beagleplink = "/app/references38/beagle_chr/plink.CHROMNAME.GRCh38.map"
    //      beagleref = "/app/references38/beagle/chrCHROMNAME.1kg.phase3.v5a_GRCh38nounref.vcf.gz"
    //    ,
    //    beagleplink = "/app/references38/beagle/plink.chrCHROMNAME.GRCh38.map"
    //    ,
  }

  def setPostFixFile(): Unit = {
    beaglejar = "/app/references38/beagle/beagle.08Feb22.fa4.jar"
    beagleref = "/app/references38/beagle/chrCHROMNAME.1kg.phase3.v5a_GRCh38nounref.vcf.gz"
    beagleplink = "/app/references38/beagle/plink.chrCHROMNAME.GRCh38.map"
    g1000prefix = "/app/references38/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chr"
    g1000alleleprefix = "/app/references38/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_allele_index_chr"
    impute_file = "/app/references38/imputation/impute_info.txt"
    problemLociFile = "/app/references38/probloci/probloci.txt.gz"
  }

