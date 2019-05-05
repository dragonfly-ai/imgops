import sbtcrossproject.CrossPlugin.autoImport.crossProject

val sharedSettings = Seq(
  version in ThisBuild := "0.2",
  scalaVersion := "2.12.6",
  organization in ThisBuild := "ai.dragonfly.code",
  scalacOptions in ThisBuild ++= Seq("-feature"),
  resolvers in ThisBuild += "dragonfly.ai" at "http://code.dragonfly.ai:8080/",
  libraryDependencies ++= Seq( "ai.dragonfly.code" %%% "img" % "0.2" ),
//  mainClass := Some("ai.dragonfly.img.ops.TestImgOps"),
  publishTo in ThisBuild := Some( Resolver.file ( "file",  new File( "/var/www/maven" ) ) )
)

val imgops = crossProject(JSPlatform, JVMPlatform)
  .settings(sharedSettings)
//  .jsSettings(scalaJSUseMainModuleInitializer := true)
