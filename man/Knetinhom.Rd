\name{Knetinhom}
\alias{Knetinhom}
\title{
  Geometrically-Corrected Inhomogeneous K Function on Network
}
\description{
  Compute the geometrically-corrected inhomogeneous \eqn{K} function
  for a point pattern on a linear network.
}
\usage{
Knetinhom(X, lambda, r = NULL, freq, ..., verbose=FALSE)
}
\arguments{
  \item{X}{
    Point pattern on a linear network (object of class \code{"lpp"}).
  }
  \item{lambda}{
    Fitted intensity of the point process. Either a numeric vector
    giving values of the fitted intensity at each data point of
    \code{X}, or an object of class
    \code{"linim"}, \code{"linfun"} or \code{"lppm"} from which the
    fitted intensity at each data point can be extracted.
  }
  \item{r}{
    Optional. Numeric vector of values of the function argument
    \eqn{r}. There is a sensible default.
  }
  \item{freq}{
    Vector of frequencies corresponding to the point events on the
    network. The length of this vector should be equal to the
    number of points on the network. The default 
    frequency is one for every point on the network.
  }
  \item{\dots}{
    Ignored.
  }
  \item{verbose}{
    Logical value indicating whether to print progress reports
    during the computation.
  }
}
\details{
  This command computes the inhomogeneous version of the
  geometrically-corrected \eqn{K} function, proposed by Ang et al
  (2012), from point pattern data on a linear network.

  The algorithm used in this computation is described in Rakshit et al (2019).

  The \pkg{spatstat} function
  \code{\link[spatstat]{linearKinhom}}
  is usable (and slightly faster) for this purpose
  for small datasets, but requires substantial amounts of memory.
  For larger datasets,
  the function \code{\link{Knetinhom}} is much more efficient.
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
   fit <- lppm(UC ~ x+y)
   r <- seq(0, 1000, length = 41)
   K <- Knetinhom(UC, lambda=fit, r = r)
}
\keyword{spatial}
\keyword{nonparametric}
