\name{decoMS2-package}
\alias{decoMS2-package}
\alias{decoMS2}
\docType{package}
\title{
Deconvolution of heterogeneous MS2 spectra collected for the purpose of identifying metabolites in untargeted metabolomics experiments
}
\description{
This package provides two functions. (1) It extracts a scan x mass matrix from MS2 scans collected with LC-MS/MS experiments in centroid mode. (2) It deconvolves contaminating fragments from the mass of interest.
}
\details{
\tabular{ll}{
Package: \tab decoMS2\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-02-24\cr
License: \tab GPL-2\cr
}
Extract the MS2 matrix using the getMS2Mat function. The most important step for successful deconvolution is to correctly specify the retention time range for getMS2Mat, by flanking the peak of interest by a few scans, without diluting out its signal by too much.

Deconvolve scans using the deconvolveSources function.
}
\author{
Igor Nikolskiy

Maintainer: Igor Nikolskiy <igor.nikolskiy@wustl.edu>
~~ The author and/or maintainer of the package ~~
}
\keyword{ package }
\examples{
# ~~ simple examples of the most important functions ~~

data(xraw)
ms2_mat= getMS2Mat(xraw, rtrange= c(27, 50))
prec= unique(xraw@msnPrecursorMz)
sources= deconvolveSources(ms2_mat, windowOffSet=c(-1, 10), cole= 20, intensityFilter= 0.01, precursors= prec)

}
