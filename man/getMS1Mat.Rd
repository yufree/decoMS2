\name{getMS1Mat}
\alias{getMS1Mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Extract dected masses from a subset of MS1 scans in xcmsRaw object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Identify masses using kernel density of centroided data.
}
\usage{
getMS1Mat(xraw, rtrange = NULL, ppm = 25, minSignal = 400)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xraw}{
%%     ~~Describe \code{xraw} here~~
an xcms raw object
}
  \item{rtrange}{
%%     ~~Describe \code{rtrange} here~~
range of retention times of scans
}
  \item{ppm}{
%%     ~~Describe \code{ppm} here~~
mass deviation between scans for obtaining EIC for each mass
}
  \item{minSignal}{
%%     ~~Describe \code{minSignal} here~~
lowest signal intensity that could be considered a mass
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{ms1_mat}{
  matrix where the rows are scans (ordered by retention time) and columns are masses
  }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Igor Nikolskiy
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (xraw, rtrange = NULL, ppm = 25, minSignal = 400)
{
    if (is.null(rtrange)) {
        scans = 1:length(xraw@scantime)
    }
    else {
        scans = which(xraw@scantime >= rtrange[1] & xraw@scantime <=
            rtrange[2])
    }
    combinedScans = do.call(rbind, sapply(scans, getScan, object = xraw))
    cons_masses = getConsensusMasses(combinedScans, ppm = ppm)
    ms1_mat = c()
    for (s in scans) {
        scan = getScan(xraw, s)
        intensities = sapply(cons_masses, function(m) {
            matches = ppm_diff(m, scan[, "mz"]) < ppm
            if (any(matches))
                return(max(scan[matches, "intensity"]))
            else return(0)
        })
        ms1_mat = rbind(ms1_mat, intensities)
    }
    colnames(ms1_mat) = cons_masses
    rownames(ms1_mat) = round(xraw@scantime[scans], 2)
    ms1_mat
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
