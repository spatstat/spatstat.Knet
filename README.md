# spatstat.Knet

### Extension to spatstat package for large data sets on a linear network

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spatstat.Knet)](http://CRAN.R-project.org/package=spatstat.Knet)
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.Knet)](https://github.com/spatstat/spatstat.Knet)

This is an _extension_ of the `spatstat` package. 

It computes the K function, pair correlation function
and inhomogeneous versions of these functions,
for point pattern on a linear network, using the
algorithms described in
S. Rakshit, A Baddeley and G. Nair (2019)
Efficient code for second order analysis of events on a linear network.
_Journal of Statistical Software_ **90** (1) 1-37.

This GitHub repository holds the *development version* of
`spatstat.Knet`. The development version is newer than the *official release*
of `spatstat.Knet` on CRAN. 

## Installing the official release

For the most recent **official release** of 
`spatstat.Knet`, see the [CRAN page](https://CRAN.R-project.org/package=spatstat.Knet). To install it, just type

```R
install.packages('spatstat.Knet')
```

## Installing the development version

The easiest way to install the **development version** of `spatstat.Knet` 
from github is through the `remotes` package:

```R
require(remotes)
install_github('spatstat/spatstat')
install_github('spatstat/spatstat.Knet')
```

If you don't have `remotes` installed you should first run

```R
install.packages('remotes')
```


