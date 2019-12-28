# @author Paul Bailey
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("DE v", utils::packageDescription("DE")$Version, "\n"))
}

globalVariables(c("dopari"))
