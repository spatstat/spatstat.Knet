#'
#' spatstat.Knet/R/First.R
#'

.onLoad <- function(...) { }

.onAttach <- function(libname, pkgname) {
  v <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                fields="Version")
  msg <- paste("spatstat.Knet", v)
  packageStartupMessage(msg)
  invisible(NULL)
}


