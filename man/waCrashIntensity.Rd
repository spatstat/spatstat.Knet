\name{waCrashIntensity}
\alias{waCrashIntensity}
\docType{data}
\title{
  A fixed-bandwidth kernel estimate of the intensity function
  for the point events in the dataset \code{wacrashes}.
}
\description{
  This is a fixed-bandwidth kernel estimate of the spatially-varying intensity
  of road traffic accidents on the road network of the state of Western
  Australia. This estimate was computed from the dataset
  \code{\link{wacrashes}} using a very fast kernel
  smoothing method proposed by Rakshit et al (2019b).
  This dataset \code{waCrashIntensity} was published and used for computing
  the inhomogeneous \eqn{K} function by Rakshit et al (2019a).
}
\usage{data(waCrashIntensity)}
\format{
  A pixel image on a linear network (object of class \code{"linim"})
  representing the intensity estimate at every location on the
  Western Australian road network.
}
\source{
  Computed by Suman Rakshit from data provided by Main Roads, Western Australia.
  Made available as part of the Western Australian
  Whole of Government Open Data Policy.
}
\references{
  Rakshit, S., Baddeley, A. and Nair, G. (2019a)
  Efficient code for second order analysis of events on a linear network.
  \emph{Journal of Statistical Software}, to appear.
  
  Rakshit, S., Davies, T., Moradi, M., McSwiggan, G., Nair, G.,
  Mateu, G. and Baddeley, A. (2019b) 
  Fast kernel smoothing of point patterns on a network using 2D convolution.
  \emph{International Statistical Review}, in press.
}
\examples{
data(waCrashIntensity)
plot(waCrashIntensity, log = TRUE)
}
\keyword{spatial}
\keyword{datasets}
\keyword{smooth}
