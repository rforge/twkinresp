#---------- process knowledge: the model relating microbial parameters to observations
kinrespModel1 <- function(
	### Respiration for given microbial parameters at given time. 
	x0		##<< initial microbial biomass (numeric scalar)
	,r0		##<< initial microbial activity (numeric scalar)
	,mumax	##<< maximum growth rate (numeric scalar)
	,time	##<< time (numeric vector)
	,lambda=0.9	##<< Ratio of growth associated (coupled) specific respiration to total specific respiration. Usually 0.9.
	,YCO2=1.5	##<< Ratio of assimilated carbon per respired carbon.  Usually 1.5.
){
	resp <- x0*(1-r0)*(1/lambda-1)*mumax/YCO2 + x0*r0*1/lambda*mumax/YCO2 * exp(mumax*time)
	### respiration at given time points (numeric vector)
}
attr(kinrespModel1,"ex") <- function(){
	data(respWutzler10)
	ds <- subset(respWutzler10, replicate==1 & experiment ==  1 & suite == 'Face' & time < 25)
	resp0 <- kinrespModel1( x0=118, r0=0.029, mumax=0.21, time=ds$time)
	plot( resp~time, data=ds)
	lines(resp0~ds$time)
}

modelPar3 <- function(
	### Translating parameter vector to arguments of Model.
	parOpt		##<< parameters: named numeric vector (x0,r0,mumax) 
	,fModel		##<< model of the form \code{function(x0,r0,mumax,...)}
	,parDefault=parOpt	
	### if parOpt contains only a subset of (x0,r0,mumax),
	### e.g. when optimizing only a subset
	### parDefault can be used to specify the other values
	,...		##<< further arguments to fModel
){
	pars <- structure( as.numeric(parDefault), names=c("x0","r0","mumax"))
	pars[names(parOpt)] <- parOpt
	fModel( pars[1], pars[2], pars[3], ...)
	### result of fModel
}

