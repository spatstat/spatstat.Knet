\name{spatstat.Knet-package}
\alias{spatstat.Knet-package}
\docType{package}
\title{
\packageTitle{spatstat.Knet}
}
\description{
\packageDescription{spatstat.Knet}
}
\details{
  This is an extension to the \pkg{spatstat} package
  for the analysis of large data sets on linear networks.

  Its main functionality is a memory-efficient algorithm
  for computing the estimate of the \eqn{K} function
  on a linear network, described in Rakshit et al (2019).

  The main functions are \code{\link{Knet}} and \code{\link{Knetinhom}}.
  These are counterparts of the functions
  \code{\link[spatstat]{linearK}}
  and 
  \code{\link[spatstat]{linearKinhom}}
  in the \pkg{spatstat} package.

  The \pkg{spatstat} functions
  \code{\link[spatstat]{linearK}}
  and 
  \code{\link[spatstat]{linearKinhom}}
  are usable (and slightly faster)
  for small datasets, but require substantial amounts of memory.
  For larger datasets,
  the functions  \code{\link{Knet}} and \code{\link{Knetinhom}}
  are much more efficient.
  
The DESCRIPTION file:
\packageDESCRIPTION{spatstat.Knet}
\packageIndices{spatstat.Knet}
}
\author{
\packageAuthor{spatstat.Knet}

Maintainer: \packageMaintainer{spatstat.Knet}
}
\references{
  Rakshit, S., Baddeley, A. and Nair, G. (2019)
  Efficient code for second order analysis of events on a linear network.
  \emph{Journal of Statistical Software} \bold{90} (1) 1--37.
  DOI: \code{10.18637/jss.v090.i01} 
}
\keyword{package}
