\name{Knet}
\alias{Knet}
\title{
  Geometrically-Corrected K Function on Network
}
\description{
  Compute the geometrically-corrected \eqn{K} function
  for a point pattern on a linear network.
}
\usage{
Knet(X, r = NULL, freq, ..., verbose=FALSE)
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument \eqn{r}. There is a sensible default.
  }
  \item{freq}{
    Vector of frequencies corresponding to the point events on the
    network. The length of this vector should be equal to the number
    of points on the network. The default 
    frequency is one for every point on the network.
  }
  \item{...}{
    Ignored.
  }
  \item{verbose}{
    A logical for printing iteration number corresponding
    to each point event on the network.
  }
}
\details{
  This command computes the geometrically-corrected \eqn{K} function,
  proposed by Ang et al (2012), from point pattern data on a linear
  network. 
  The algorithm used in this computation is discussed in Rakshit et al (2019).

  The \pkg{spatstat} function \code{\link[spatstat]{linearK}}
  is usable (and slightly faster) for the same purpose
  for small datasets, but requires substantial amounts of memory.
  For larger datasets, the function  \code{\link{Knet}} is much more efficient.
}
\value{
  Function value table (object of class \code{"fv"}).
}
\author{
  Suman Rakshit (modified by Adrian Baddeley)
}
\references{
  Ang, Q.W., Baddeley, A. and Nair, G. (2012)
  Geometrically corrected second-order analysis of 
  events on a linear network, with applications to
  ecology and criminology.
  \emph{Scandinavian Journal of Statistics} \bold{39}, 591--617.
  
  Rakshit, S., Baddeley, A. and Nair, G. (2019)
  Efficient code for second order analysis of events on a linear network.
  \emph{Journal of Statistical Software} \bold{90} (1) 1--37.
  DOI: \code{10.18637/jss.v090.i01} 
}
\examples{
   UC <- unmark(chicago)
   r <- seq(0, 1000, length = 41)
   K <- Knet(UC, r = r)
}
\keyword{spatial}
\keyword{nonparametric}
