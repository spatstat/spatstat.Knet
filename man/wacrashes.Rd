\name{wacrashes}
\alias{wacrashes}
\docType{data}
\title{
  Road Accidents in Western Australia
}
\description{
  This dataset gives the spatial locations of all road accidents
  recorded in the state of Western Australia for the year 2011, on the
  state road network.
  
  These data were published and analysed in Rakshit et al (2019).
}
\usage{data(wacrashes)}
\format{
  A object of class \code{"lpp"} representing the spatial point pattern
  of accident locations on the network of roads in Western Australia.
}
\details{
  The road network has 88,512 intersections and 115,169 road segments.
  The spatial coordinates are expressed in metres, and the
  total network length is 97,165,540 metres (97,165 km).
  The number of accident locations on the network is 14,562.
}
\source{
  Main Roads, Western Australia.
  Made available as part of the Western Australian
  Whole of Government Open Data Policy. 
}
\references{
  Rakshit, S., Baddeley, A. and Nair, G. (2019)
  Efficient code for second order analysis of events on a linear network.
  \emph{Journal of Statistical Software} \bold{90} (1) 1--37.
  DOI: \code{10.18637/jss.v090.i01} 
}
\examples{
  data(wacrashes)
  wacrashes
  summary(wacrashes)
  plot(wacrashes, cols="red", cex=0.5)
}
\keyword{spatial}
\keyword{datasets}
