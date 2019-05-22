#'
#' spatstat.Knet/R/Knetinhom.R
#'
#' Copyright (C) 2015-2019 Suman Raskhit and Adrian Baddeley
#' GNU Public Licence GPL >= 2
#'

Knetinhom <- function(X, lambda, r = NULL, freq, ..., verbose = FALSE) {
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
  
  if(is.numeric(lambda)) {
    check.nvector(lambda, npoints(X), things="points")
    lambdavalues <- lambda
  } else {
    ## extract lambda values
    if(inherits(lambda, "linim")) {
      lambdavalues <- lambda[X, drop=FALSE]
    } else if(inherits(lambda, "linfun")) {
      lambdavalues <- lambda(X)
    } else if(is.lppm(lambda)) {
      lambdavalues <- predict(lambda, locations=coords(X))
    } else stop("Format of argument lambda is not understood")
    ok <- check.nvector(lambdavalues, npoints(X),
                        fatal=FALSE, warn=TRUE)
      if(!ok) stop("Incorrect format for lambda values")
  }
  ## assemble data 
  ## points: x, y, seg, tp, lambda, freq
  nX <- npoints(X)
  if(missing(freq)){
    df <- cbind(coords(X), data.frame(lambda=lambdavalues, freq=1))
  }else{
    f <- as.integer(freq)
    if(length(f) == nX & all(f) > 0){
      df <- cbind(coords(X),data.frame(lambda=lambdavalues, freq=f))
    }else{
      stop("freq vector is not correctly specified.")
    }
  }
  ## tweak
  
  ## TWo cases needs to be considered in order to aggregate the "lambda"-values.
  ## Case-I Aggregate points with tp = 0 or tp = 1. In this case, "freq" should add for all identical
  ## points, but we do not add the lambda values.
  #df[which(df$tp == 1.0), "tp"] <- 1-Eps
  #df[which(df$tp == 0.0), "tp"] <- Eps
  dfFreq <- aggregate(freq ~ seg + tp, data=df, FUN="sum")
  dfLamb <- aggregate(lambda ~ seg + tp, data=df, FUN="mean")
  df <- data.frame(dfFreq,lambda=dfLamb$lambda)
  ## Case-II Add the lambda values for distinct points on the same segment with 0 < tp < Eps or 1-Eps < tp < 1.
  Eps <- 1e-4
  df$tp <- pmax(Eps, pmin(1-Eps, df$tp))
  # sort by seg and by tp within seg
  dfGrpd <- aggregate(cbind(freq,lambda) ~ seg + tp, data=df, FUN="sum")
  ord <- with(dfGrpd, order(seg, tp))
  df <- dfGrpd[ord, , drop=FALSE]
  Seg  <- df$seg
  Tp   <- df$tp
  Freq <- df$freq
  Lamb <- df$lambda
  ## linear network
  L <- as.linnet(X)
  noVert <- npoints(vertices(L))
  PL <- as.psp(L)
  noEdge <- nobjects(PL)
  EdgeLengths <- lengths.psp(PL)
  nX <- nrow(df)
  Vert1 <- L$from
  Vert2 <- L$to
  ## call C code
  z<- .C("I_createGraphNet",
         no_of_vertices=as.integer(noVert),
         no_of_edges=as.integer(noEdge),
         no_of_crashes=as.integer(nX),
         crash_seg=as.integer(Seg), 
         crash_tp=as.double(Tp),
         crash_freq=as.integer(Freq),
         crash_lambda=as.double(Lamb),
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
               ylab = quote(K[L, inhom](r)),
               yexp = quote(K[list(L, "inhom")](r)),
               valu = "est",
               labl = c("r", "{%s[%s]^{pois}}(r)", "{hat(%s)[%s]}(r)"),
               desc = c("distance argument r",
                        "theoretical Poisson %s",
                        "estimate of %s"),
               fname = c("K", "list(L, inhom)"),
               fmla = . ~ r,
               unitname = unitname(X)
  )
  return(result)
}

