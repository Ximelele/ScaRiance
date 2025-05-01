ThisBuild / version := "0.1.0-SNAPSHOT"
ThisBuild / scalaVersion := "3.3.5"


lazy val root = (project in file("."))
  .settings(
    name := "Battenberg",
  )


libraryDependencies ++= Seq(

  ("org.apache.spark" %% "spark-sql" % "3.5.5").cross(CrossVersion.for3Use2_13),
  ("org.apache.spark" %% "spark-core" % "3.5.5").cross(CrossVersion.for3Use2_13),
  ("org.scala-lang.modules" %% "scala-parallel-collections" % "1.2.0").cross(CrossVersion.for3Use2_13)

)


fork := true
run / javaOptions ++= Seq("-Xmx5G", "-XX:+UseG1GC")

