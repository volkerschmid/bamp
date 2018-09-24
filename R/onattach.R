.onAttach<-function(libname, pkgname)
{
  packageStartupMessage(paste0("BAMP ver. ", utils::packageVersion("bamp")))
}
