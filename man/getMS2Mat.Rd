\name{getMS2Mat}
\alias{getMS2Mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Extract a scan x mass matrix from xcmsRaw object(s).
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
For the defined scans, this function identifies consensus masses based on a minimal signal intensity cut-off. Then a data.frame is returned with EICs in columns for each identified mass.
}
\usage{
getMS2Mat(xraw, region = NULL, precursors = NULL, rtrange = NULL, collisionEnergy = NULL, ppm = 30, minSignal = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xraw}{
%%     ~~Describe \code{xraw} here~~
xcmsRaw object or list of xcmsRaw objects
}
  \item{region}{
%%     ~~Describe \code{region} here~~
rows from a peak table where the feature of interest was detected (optional)
}
  \item{precursors}{
%%     ~~Describe \code{precursors} here~~
the precursors that the mass spec targeted
}
  \item{rtrange}{
%%     ~~Describe \code{rtrange} here~~
range of scans to use
}
  \item{collisionEnergy}{
%%     ~~Describe \code{collisionEnergy} here~~
the collision energy of the scans to use
}
  \item{ppm}{
%%     ~~Describe \code{ppm} here~~
ppm error mass variability between scans
}
  \item{minSignal}{
%%     ~~Describe \code{minSignal} here~~
smallest intensity that could initiate a mass to include as column
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
if precursors, rtrange, or collsionEnergy are left out or NULL, all values present in xraw are used
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{ms2_mat}{
  a data.frame containing a row for each scan, with:
  xcmsRaw:  index of xcmsRaw file
  scan:         scan number
  rt:             retention time of scan
  ce:            collision energy of scan
  prec:         precursor that was targeted
  and a column for every mass.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Igor Nikolskiy
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (xraw, region = NULL, precursors = NULL, rtrange = NULL,
    collisionEnergy = NULL, ppm = 30, minSignal = 100)
{
    if (!is.null(region)) {
        rtrange = c(min(region$rtmin) - 5, max(region$rtmax) +
            5)
    }
    if (class(xraw) == "xcmsRaw") {
        xraw = c(xraw)
    }
    scans = lapply(xraw, selectScans, precursors = precursors,
        rtrange = rtrange, collisionEnergy = collisionEnergy)
    combinedScans = combineScans(xraw, scans, minSignal)
    if (nrow(combinedScans) > 2) {
        cons_masses = getConsensusMasses(combinedScans, ppm = ppm)
        ms2_mat = createMatrix(xraw, scans, cons_masses, ppm = ppm)
        return(ms2_mat)
    }
    else {
        cat("ms2 has no masses with intensity above: ", minSignal,
            "\n")
    }
  }
}
