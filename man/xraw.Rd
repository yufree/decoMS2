\name{xraw}
\alias{xraw}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Experimental data to deconvolve cysteine
}
\description{
Experimental data to deconvolve cysteine
}
\usage{data(xraw)}
\format{
  The format is:
Formal class 'xcmsRaw' [package "xcms"] with 20 slots
  ..@ env                  :<environment: 0x106258430>
  ..@ tic                  : num [1:62] 2106 1929 237638 960298 662368 ...
  ..@ scantime             : num [1:62] 2.42 2.91 5.91 8.91 11.91 ...
  ..@ scanindex            : int [1:62] 0 80 152 1319 3448 5624 7686 9899 12004 14092 ...
  ..@ polarity             : Factor w/ 3 levels "negative","positive",..: 1 1 1 1 1 1 1 1 1 1 ...
  ..@ acquisitionNum       : int [1:62] 2428 2915 5916 8919 11919 14921 17924 20925 23927 26928 ...
  ..@ profmethod           : chr "bin"
  ..@ profparam            : list()
  ..@ mzrange              : num [1:2] 30 200
  ..@ gradient             : logi[0 , 0 ]
  ..@ msnScanindex         : int [1:236] 0 71 206 314 484 793 1166 1359 1624 2149 ...
  ..@ msnAcquisitionNum    : int [1:236] 3911 4409 4908 5407 6911 7410 7909 8409 9914 10412 ...
  ..@ msnPrecursorScan     : int [1:236] NA NA NA NA NA NA NA NA NA NA ...
  ..@ msnLevel             : int [1:236] 2 2 2 2 2 2 2 2 2 2 ...
  ..@ msnRt                : num [1:236] 3.9 4.4 4.9 5.4 6.9 ...
  ..@ msnPrecursorMz       : num [1:236] 121 121 114 114 121 121 114 114 121 121 ...
  ..@ msnPrecursorIntensity: num [1:236] 2646 7695 5013 10508 25815 ...
  ..@ msnPrecursorCharge   : int [1:236] 1 1 1 1 1 1 1 1 1 1 ...
  ..@ msnCollisionEnergy   : num [1:236] 20 0 20 0 20 0 20 0 20 0 ...
  ..@ filepath             : chr "../data/aa_deconvolution/shiftingWindows/set1_cys.mzML"
}
\examples{
data(xraw)
## maybe str(xraw) ; plot(xraw) ...
}
\keyword{datasets}
