#'
#' spatstat.Knet/R/Knet.R
#'
#' Copyright (C) 2015-2019 Suman Rakshit and Adrian Baddeley
#' GNU Public Licence GPL >= 2
#'

Knet <- function(X, r = NULL, freq, ..., verbose = FALSE) {
  ## validate data
  stopifnot(is.lpp(X))
  
  if(missing(r) || is.null(r)){
    Mndist <- 0
    Mxdist <- 2* mean(nndist(as.ppp(X)))
    noGrid <- 41
  }else{
    stopifnot(is.numeric(r))
    Mndist <- min(r)
    Mxdist <- max(r)
    noGrid <- length(r)
  }
  stopifnot(noGrid >= 2)
  maxR <- Mxdist + (Mxdist/100)
  ## assemble data 
  nX <- npoints(X)
  ## points: x, y, seg, tp, freq
  if(missing(freq)){
    df <- cbind(coords(X), data.frame(freq=rep(1, nX)))
  }else{
    f <- as.integer(freq)
    if(length(f) == nX & all(f) > 0){
      df <- cbind(coords(X),data.frame(freq=f))
    }else{
      stop("freq vector is not correctly specified.")
    }
  }
  ## tweak
  Eps <- 1e-4
  df$tp <- pmax(Eps, pmin(1-Eps, df$tp))
  # sort by seg and by tp within seg
  dfGrpd <- aggregate(freq ~ seg + tp, data=df, FUN="sum")
  ord <- with(dfGrpd, order(seg, tp))
  df <- dfGrpd[ord, , drop=FALSE]
  Seg  <- df$seg
  Tp   <- df$tp
  Freq <- df$freq
  nX <- nrow(df)
  ## linear network
  L <- as.linnet(X)
  noVert <- npoints(vertices(L))
  PL <- as.psp(L)
  noEdge <- nobjects(PL)
  EdgeLengths <- lengths.psp(PL)
  Vert1 <- L$from
  Vert2 <- L$to
  ## call C code
  z<- .C("createGraphNet",
         no_of_vertices=as.integer(noVert),
         no_of_edges=as.integer(noEdge),
         no_of_crashes=as.integer(nX),
         crash_seg=as.integer(Seg), 
         crash_tp=as.double(Tp),
         crash_freq=as.integer(Freq),
         vert_id1=as.integer(Vert1),
         vert_id2=as.integer(Vert2), 
         edge_length=as.double(EdgeLengths),
         MAX_Distance = as.double(Mxdist),
         MIN_Distance = as.double(Mndist),
         no_of_distance = as.integer(noGrid),
         max_r = as.double(maxR),
         verboseIter = as.integer(verbose),
         kvalue = as.double(numeric(noGrid)),
         PACKAGE="spatstat.Knet")
  
  r <- seq(0, maxR, length.out = noGrid)
  df <- data.frame(r=r, theo=r, est=z$kvalue)
  result <- fv(df,
               argu = "r",
               ylab = quote(K[L](r)),
               valu = "est",
               labl = c("r", "{%s[%s]^{pois}}(r)", "{hat(%s)[%s]}(r)"),
               desc = c("distance argument r",
                        "theoretical Poisson %s",
                        "estimate of %s"),
               fname = c("K", "L"),
               fmla = . ~ r,
               unitname = unitname(X)
               )
  return(result)
}

