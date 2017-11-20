.onAttach<-function(libname, pkgname)
{
  packageStartupMessage(paste0("bamp ver. ", utils::packageVersion("bamp")))
}
