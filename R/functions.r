


library(RColorBrewer)
library(xcms)
library(gplots)
library(nnls)



##############################
######## functions to get ms2_mat ###
##############################

selectPeaks <- function(peaks, rtrange, mzrange){
	peaks= as.data.frame(peaks)
	p= peaks$rt>rtrange[1] &
		peaks$rt<rtrange[2] &
		peaks$mz>mzrange[1] &
		peaks$mz<mzrange[2]
	peaks[p,]
}


getMS2Scan <- function(object, scan=NULL, acqnum=NULL, normalize = FALSE, filtersmall=TRUE) {
	if (is.null(scan))
		scan <- which(object@msnAcquisitionNum == acqnum)
	## handle last spectrum
	if (scan == length(object@msnScanindex)) {
		followingScanIndex <- length(object@env$msnMz)
	} else {
		followingScanIndex <- object@msnScanindex[scan+1]
	}
	## handle empty spectra
	if (object@msnScanindex[scan] == length(object@env$msnMz) ||
		object@msnScanindex[scan] == followingScanIndex) {
		warning("empty scan")
		return()
	}

	idx <- (object@msnScanindex[scan]+1):min(followingScanIndex,length(object@env$msnMz), na.rm=TRUE)

	points <- cbind(object@env$msnMz[idx], object@env$msnIntensity[idx])
	colnames(points) <- c("mz","intensity")

	if (normalize)
		points[,"intensity"] <- (points[,"intensity"] / max(points[,"intensity"])) * 100

	if (filtersmall)
		points <- points[points[,"intensity"] > 0.001,,drop=FALSE]

	points
}


selectScans <- function(xraw, precursors= NULL, rtrange= NULL,
                        collisionEnergy= NULL){
	if (is.null(collisionEnergy)){
		collisionEnergy= unique(xraw@msnCollisionEnergy)
	}
	if (is.null(precursors)){
		precursors= unique(xraw@msnPrecursorMz)
	}
	if (is.null(rtrange)){
		scans= which( xraw@msnPrecursorMz %in% precursors &
		xraw@msnCollisionEnergy %in% collisionEnergy)
	}else{
		scans= which( xraw@msnRt >= rtrange[1] &
		xraw@msnRt <= rtrange[2] &
		xraw@msnPrecursorMz %in% precursors &
		xraw@msnCollisionEnergy %in% collisionEnergy)
	}
	scans
}


combineScans <- function(xraw, scans, minSignal= 100){
	pl= c()
	for (x in 1: length(xraw)){
		for (i in 1:length(scans[[x]])){
			scan= getMS2Scan(xraw[[x]], scans[[x]][i])
			intenseEnough= (scan[,'intensity'] >= minSignal)
			pl= rbind(pl, scan[intenseEnough, , drop= F])
		}
	}
	pl
}


filterRipple <- function(d, maxRange= 0.1){
	masses= c()
	remainingMasses= 1:nrow(d)
	while( any(remainingMasses) ){
		m= d[remainingMasses[ 1 ], 1]
		matches= which(abs(m-d[remainingMasses,1])<=maxRange)
		masses= rbind(masses, c(mz= m, intensity= d[remainingMasses[ 1 ], 2]))
		remainingMasses= remainingMasses[-matches]
	}
	masses
}


findPeakMasses <- function(kd, maxRange= 0.105, minInt= 0.01, round_to= 4){
	d= cbind(kd$x, kd$y)
	d= d[order(d[,2], decreasing=T),]
	d= d[d[,2]>=d[1,2]*minInt, ]
	masses= filterRipple(d, maxRange-.04)
	masses= filterRipple(masses, maxRange+0.01)
	masses= sort(round(masses[,1], round_to))
	masses
}


getConsensusMasses <- function(combinedScans, ppm= 30, round= 4, deisotope= FALSE){
	pl= combinedScans
	n= diff(range(pl[,1]))*1000
	kd= density(pl[,1], weights= pl[,2]/sum(pl[,2]), bw= 0.001, n= n)
	masses= findPeakMasses(kd)
	# if (deisotope) {
	# 	isotopes= which(abs(diff(masses)-1.00335)/masses[1:(length(masses)-1)]*1e6 < ppm) + 1
	# 	masses= masses[-isotopes]
	# }
	masses
}


ppm_diff <- function(mass, masses){
	abs(mass-masses)/mass*1e6
}


createMatrix <- function(xraw, scans, cons_masses, ppm= 30){
	mat= c()
	for (x in 1:length(xraw)){
		# cat(x, '\n')
		for (s in scans[[x]]){
			ms2_scan= getMS2Scan(xraw[[x]], s)
			intensities= sapply(cons_masses, function(m){
													matches= ppm_diff(m, ms2_scan[,'mz'])<ppm
													if (any(matches)) return(max(ms2_scan[matches, 'intensity']))
													else return(0)
												})
			names(intensities)= cons_masses
			row= c(	xcmsRaw= x,
						scan= s,
						rt= xraw[[x]]@msnRt[s],
						ce= xraw[[x]]@msnCollisionEnergy[s],
						prec= xraw[[x]]@msnPrecursorMz[s],
						intensities
						)
			mat= rbind(mat, row)
		}
	}
	rownames(mat)= NULL
	as.data.frame(mat)
}


getMS2Mat <- function(	xraw,
							region= NULL,
							precursors= NULL,
							rtrange= NULL,
							collisionEnergy= NULL,
							ppm= 30,
							minSignal= 100){
	# '''
	# This function accepts a list of xraw files and parameters specifying what portion of those files to work with.
	# It returns a matrix with masses as columns and scans as rows.
	# '''
	if (! is.null(region) ) {
		rtrange= c( min(region$rtmin)-5, max(region$rtmax)+5 )
		## possibly choose from the precursors here, if anyone cares
	}
	if (class(xraw)=='xcmsRaw'){
		xraw= c(xraw)
	}
	scans= lapply(	xraw,
					selectScans,
					precursors= precursors,
					rtrange= rtrange,
					collisionEnergy= collisionEnergy)
	combinedScans= combineScans(xraw, scans, minSignal)
	if (nrow(combinedScans)>2){
		cons_masses= getConsensusMasses(	combinedScans,
												ppm= ppm)
		ms2_mat= createMatrix(xraw, scans, cons_masses, ppm= ppm)
		return(ms2_mat)
	}else{
		cat('ms2 has no masses with intensity above: ', minSignal, '\n')
	}
}


getMS1Mat <- function(xraw, rtrange= NULL, ppm= 25, minSignal= 400){
	if (is.null(rtrange)){
 		scans= 1:length(xraw@scantime)
	}else{
		scans= which(xraw@scantime >= rtrange[1] & xraw@scantime <= rtrange[2])
	}
	combinedScans= do.call(rbind, sapply(scans, getScan, object= xraw))
	# combinedScans= combinedScans[combinedScans[,2] >= minSignal, ]
	cons_masses= getConsensusMasses(	combinedScans,
											ppm= ppm)

	ms1_mat= c()
	for (s in scans){
		scan= getScan(xraw, s)
		intensities= sapply(cons_masses, function(m){
												matches= ppm_diff(m, scan[,'mz'])<ppm
												if (any(matches)) return(max(scan[matches, 'intensity']))
												else return(0)
											})
		ms1_mat= rbind(ms1_mat, intensities)
	}
	colnames(ms1_mat)= cons_masses
	rownames(ms1_mat)= round(xraw@scantime[scans], 2)
	ms1_mat
}


########################################
####### Functions to deconvolve sources     #######
########################################
subset_ms2_mat <- function(ms2, cole= NULL, prec= NULL, xr=NULL, minColPosPortion= 0){
	if (is.null(cole)) cole= unique(ms2$ce)
	if (is.null(prec)) prec= unique(ms2$prec)
	if (is.null(xr)) xr= unique(ms2$xcmsRaw)

	rows= ms2$prec %in% prec & ms2$ce==cole & ms2$xcmsRaw %in% xr
	s= ms2[rows, grep("[0-9]+.[0-9]+", colnames(ms2)), drop= F]
	rownames(s)= round(ms2$rt[rows], 3)
	cols= apply(s, 2, function(x) sum(x>100)/nrow(s)>=minColPosPortion)

	s[,cols]
}

findIsotopeMasses <- function(masses, ppm= 30){
	n= length(masses)
	isotopes= c()
	for (i in 1:n){
		for (j in i:n){
			d= abs(masses[j]-masses[i]-1.003355)/masses[j]*1e6<ppm
			if (d){
				# cat(d, i, j, '\n')
				isotopes= c(isotopes, j)
			}
		}
	}
	isotopes
}

interpolatePrecursor <- function(intensity, observedRts, targetRts, smoothFunc= 'spline'){
	observedRts= as.double(observedRts)
	targetRts= as.double(targetRts)
	if (smoothFunc=='spline'){
		f= smooth.spline(intensity ~ observedRts)
		return(predict(f, targetRts)$y)
	}else if (smoothFunc=='loess'){
		f= loess(intensity ~ observedRts, span= 0.15)
		return(predict(f, targetRts))
	}else if (smoothFunc=='linear'){
		return(approx(observedRts, intensity, targetRts)$y)
	}
}

deconvolveSources <- function(ms2_mat, precursors, windowOffSet= c(-1, 9), deisotope= T, intensityFilter= 0.05, cole= 20){
	precMat= subset_ms2_mat(ms2_mat, prec= precursors, cole= 0)
	cols= grep("[0-9]+.[0-9]+", colnames(precMat))
	masses= as.numeric(colnames(precMat)[cols])
	possiblePrecs= min(precursors) + windowOffSet[1] < masses & max(precursors)+windowOffSet[2]>masses
	cols= cols[possiblePrecs]
	masses= as.numeric(colnames(precMat)[cols])
	if (deisotope) {
		isos= findIsotopeMasses(masses)
		if (any(isos)) cols= cols[-isos]
	}
	interpolated= c()
	observations= c()
	for(xr in unique(ms2_mat$xcmsRaw)){
		for(p in unique(ms2_mat$prec)){
			om= subset_ms2_mat(ms2_mat, prec= p, cole= cole, xr= xr)
			pm= subset_ms2_mat(ms2_mat, prec= p, cole= 0, xr= xr)
			loc= apply(pm[,cols, drop= F], 2, function(x){
						interpolatePrecursor(x, rownames(pm), rownames(om))
						})
			interpolated= rbind(interpolated, loc)
			observations= rbind(observations, om)
		}
	}
	cols= max(interpolated)*intensityFilter < apply(interpolated, 2, max)
	interpolated= interpolated[,cols, drop= F]

	# sources= ginv(interpolated) %*% observations
	sources= sapply(apply(observations, 2, nnls, A= interpolated), function(x) x$x)
	if(! is.null(dim(sources))){
		rownames(sources)= colnames(interpolated)
	}
	sources
}



########################################
####### plotting functions and etc. #############
########################################

filterFragments <- function(ms2_mat){
	mzs= names(ms2_mat)
	mzs_int= floor(as.double(mzs))
	toKeep=c()
	for (m in unique(mzs_int)){
		matches= mzs[ mzs_int == m ]
		if (length(matches)==1){
			toKeep= c(toKeep, matches)
		}else{
			cs= ms2_mat[matches]
			toKeep= c(toKeep, matches[ cs == max(cs) ])
		}
	}
	toReturn= ms2_mat[toKeep]
	toReturn
}

getTopMasses <- function(sr, N=20, filter= TRUE){
	if (filter) sr= filterFragments(sr)
	masses= as.double(names(sr))
	if (length(sr)<=N){
		top= 1:length(sr)
	}else{
		top= which(sr>sort(sr)[length(sr)-N])
	}
	top= sr[top]/max(sr)
	top= cbind(as.double(names(top)), top)
	colnames(top)= c('mz', 'relative intensity')
	rownames(top)= NULL
	top
}

plot_ms2_mat <- function(ms2, title= ''){
	ms2= as.matrix(ms2)
	heatmap.2(ms2,
				trace= 'none',
				dendrogram= 'none',
				Rowv= NA,
				Colv= NA,
				scale= 'none',
				cexCol= .75,
				cexRow= .75,
				col= colorRampPalette(brewer.pal(9, 'BuPu'))(100),
				main= title)
}

plotEICs <- function(mat, path= NULL){
	if (! is.null(path)) pdf(path)
	for (i in 1:ncol(mat)){
		plot( mat[,i], type= 'b', main= colnames(mat)[i],
		ylab= 'intensity', xlab= 'scan')
	}
	if (! is.null(path)) dev.off()
}

plot_source <- function(sources, i){
	plot(getTopMasses(sources[i,]), type= 'h', main= rownames(sources)[i])
}








