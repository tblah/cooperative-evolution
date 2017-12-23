import Dependencies._

lazy val root = (project in file(".")).
  settings(
    inThisBuild(List(
      organization := "com.example",
      scalaVersion := "2.12.3",
      version      := "0.1.0-SNAPSHOT"
    )),
    name := "cooperative-evolution",
    libraryDependencies  ++= Seq(
    // Last stable release
    "org.scalanlp" %% "breeze" % "0.13.2",
    
    // Native libraries are not included by default. add this if you want them (as of 0.7)
    // Native libraries greatly improve performance, but increase jar sizes. 
    // It also packages various blas implementations, which have licenses that may or may not
    // be compatible with the Apache License. No GPL code, as best I know.
    "org.scalanlp" %% "breeze-natives" % "0.13.2",
    
    // The visualization library is distributed separately as well.
    // It depends on LGPL code
    "org.scalanlp" %% "breeze-viz" % "0.13.2"
  ),
  resolvers += "Sonatype Releases" at "https://oss.sonatype.org/content/repositories/releases/",
  libraryDependencies += "org.sameersingh.scalaplot" % "scalaplot" % "0.0.4",
  // additional heap
  javaOptions += "-Xmx3G",
  fork := true,
  )
