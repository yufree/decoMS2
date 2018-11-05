\name{deconvolveSources}
\alias{deconvolveSources}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Deconvolves the mixed scans in ms2_mat returned by getMS2Mat
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Uses non-negative least squares to find the most likely portion of each observed mass that was contributed by each precursor in each scan
}
\usage{
deconvolveSources(ms2_mat, precursors, windowOffSet = c(-1, 9), deisotope = T, intensityFilter = 0.05, cole = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms2_mat}{
%%     ~~Describe \code{ms2_mat} here~~
data.frame returned by getMS2Mat
}
  \item{precursors}{
%%     ~~Describe \code{precursors} here~~
the precursors that were targeted
}
  \item{windowOffSet}{
%%     ~~Describe \code{windowOffSet} here~~
window offset for the instrument that gives the bounds of the isolation window. Used to define which columns could potentially be precursors.
}
  \item{deisotope}{
%%     ~~Describe \code{deisotope} here~~
remove potential isotope precursors
}
  \item{intensityFilter}{
%%     ~~Describe \code{intensityFilter} here~~
minimal relative intensity for a precursor
}
  \item{cole}{
%%     ~~Describe \code{cole} here~~
collsiion energy that contains non0V scans
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Most useful parameter to tinker with is the intensityFilter
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{sources}{
    a matrix with a row for each source and column for each mass. Plot each row with plot_source.
}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
function (ms2_mat, precursors, windowOffSet = c(-1, 9), deisotope = T,
    intensityFilter = 0.05, cole = 20)
{
    precMat = subset_ms2_mat(ms2_mat, prec = precursors, cole = 0)
    cols = grep("[0-9]+.[0-9]+", colnames(precMat))
    masses = as.numeric(colnames(precMat)[cols])
    possiblePrecs = min(precursors) + windowOffSet[1] < masses &
        max(precursors) + windowOffSet[2] > masses
    cols = cols[possiblePrecs]
    masses = as.numeric(colnames(precMat)[cols])
    if (deisotope) {
        isos = findIsotopeMasses(masses)
        if (any(isos))
            cols = cols[-isos]
    }
    interpolated = c()
    observations = c()
    for (xr in unique(ms2_mat$xcmsRaw)) {
        for (p in unique(ms2_mat$prec)) {
            om = subset_ms2_mat(ms2_mat, prec = p, cole = cole,
                xr = xr)
            pm = subset_ms2_mat(ms2_mat, prec = p, cole = 0,
                xr = xr)
            loc = apply(pm[, cols, drop = F], 2, function(x) {
                interpolatePrecursor(x, rownames(pm), rownames(om))
            })
            interpolated = rbind(interpolated, loc)
            observations = rbind(observations, om)
        }
    }
    cols = max(interpolated) * intensityFilter < apply(interpolated,
        2, max)
    interpolated = interpolated[, cols, drop = F]
    sources = sapply(apply(observations, 2, nnls, A = interpolated),
        function(x) x$x)
    if (!is.null(dim(sources))) {
        rownames(sources) = colnames(interpolated)
    }
    sources
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
