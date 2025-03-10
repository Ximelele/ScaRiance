import org.apache.spark.sql.SparkSession

@main
def main(): Unit =
  val battenberg = Battenberg("/data/I002.007.WGS.C.bam", "/data/I002.007.WGS.T.bam")

  //  val sparkSession = SparkSession.builder()
  //    .appName("Battenberg")
  //    .master("local[*]")
  //    .getOrCreate()

  battenberg.run()
