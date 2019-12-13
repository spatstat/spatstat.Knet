#' tests/alltests.R
#' Tests for the spatstat.Knet package

require(spatstat.Knet)
local({
  #' network data not ordered (from < to)
  L <- simplenet
  a <- L$from[c(FALSE,TRUE)]
  L$from[c(FALSE,TRUE)] <- L$to[c(FALSE,TRUE)]
  L$to[c(FALSE,TRUE)] <- a
  X <- runiflpp(20, L)
  K <- Knet(X)
  if(all(K$est == 0)) stop("Knet failed, when network data are not in order")
})
