fitKinrespBetaReplicate <- function(
	### Fitting the beta-model with estimating the log of beta0, so that exp will be positive.
	x
	,y
	, weights=NULL
){
	##details<< \describe{\item{Functions related to fitting the beta-form of the model.}{
	## \describe{
	## \item{\code{\link{modelKinrespBeta}} }{ Model form with simple coefficients beta0..beta2.  }
	## \item{\code{\link{calcKinrespCoef}} }{ Calulate Microbial kinetic parameters from given beta coefficients.  }
	## \item{\code{\link{coefKinrespBeta.default}} }{ Transform coefficients of Beta-form model from transformed scale to original scale.  }
	## }
	##}}
	
	##seealso<< 
	## \code{\link{twKinresp}}
	## , \code{\link{kinrespGrowthphaseReplicate}}
	tmp.start <- coef(lm(log(y-0.99*min(y)) ~ x))
	# fit the log-Transformed beta01 and beta21
	tmp.fit <- try(gnls( y ~ exp(beta0l) + beta1 * exp( exp(beta2l) * x)
			, params=beta0l+beta1+beta2l~1
			, start=list(beta0l=log(min(y)),beta1=exp(tmp.start[1]),beta2l=log(tmp.start[2]))
			, weights=weights
			, data=data.frame(x=x, y=y)
		))
	if(inherits(tmp.fit,"try-error")){
		# fit with fixed activity ratio of 1
		tmp.m <- .fitExpModelConstr(x,y,weights)
		#tmp.coef <- coefKinrespBeta(coef(tmp.m)) # tmp.cf <- confint(tmp.m); confintKinrespBeta(tmp.cf) 
		#lines(fitted(tmp.m)~x)
		tmp.start <- .coefKinrespBetaLogStart(coef(tmp.m))
		if( is.infinite(tmp.start[1])) tmp.start[1] <- -1e8
		tmp.fit <- try(gnls( y ~ exp(beta0l) + beta1 * exp( exp(beta2l) * x)
				, params=beta0l+beta1+beta2l~1
				, start=tmp.start
				, weights=weights
				, data=data.frame(x=x, y=y)
			))
		if(inherits(tmp.fit,"try-error"))
			tmp.fit <- tmp.m
	}
	tmp.fit
}

modelKinrespBeta <- function(
	### Calculate respiration for time x based on coefficients beta_i.
	x		##<< time
	,betai	##<< named numeric vector (beta0,beta1,beta2)
){
	##seealso<< 
	## \code{\link{fitKinrespBetaReplicate}}
	betai["beta0"] + betai["beta1"]*exp(betai["beta2"]*x) 
	### \code{beta0 + beta1*exp(beta2*x)}
}


.fitExpModel <- function(
	### fit y ~ beta0 + beta1 * exp( beta2 * x)
	x
	,y
	,weights=NULL
){
	# fitExpModel
	tmp.start <- coef(lm(log(y-0.99*min(y)) ~ x))
	# coefficients at original scale, might be negative
	tmp.fit <-gnls( y ~ beta0 + beta1 * exp( beta2 * x)
		, params=beta0+beta1+beta2~1
		, start=list(beta0=min(y),beta1=exp(tmp.start[1]),beta2=tmp.start[2])
		, weights=weights
		, data=data.frame(x=x, y=y)
	)
	tmp.fit
	### result of \code{\link{gnls}}
}

.fitExpModelConstr <- function(
	### Invokes fitExpModel to apply the beta_0..beta_2 form to the x,y, data.
	x
	,y
	,weights=NULL
	, r0=NULL
){
	# fitExpModelConstr
	#
	##details<<
	## if beta0 is estimated <0 then fit beta1..2 form (beta0=0) to the xy data
	## y ~ beta1 * exp( beta2 * x)
	if( is.null(r0) ){
		tmp.fit <- .fitExpModel(x,y,weights)
		if( coef(tmp.fit)["beta0"] < 0) r0 = 0
	}
	if( !is.null(r0) ) {
		# hard-set beta0 to zero then 
		tmp.start <- coef(lm(log(y) ~ x))
		tmp.fit <-gnls( y ~ beta1 * exp( beta2 * x)
			, params=beta1+beta2~1
			, start=list(beta1=exp(tmp.start[1]),beta2=tmp.start[2])
			, weights=weights
			, data=data.frame(x=x, y=y)
		)
	}
	tmp.fit
}

.coefKinrespBetaLogStart <- function(
	### Transform beta to log scale.
	tmp.coef	## beta coefficients at original scale
){
	# coefKinrespBetaLogStart
	##details<< 
	## Replaces the first coefficient ("beta0") by its logarithm and renames to "beta0l".
	## Replaces the third coefficient ("beta2") by its logarithm and renames to "beta2l".
	tmp.coef <- coefKinrespBeta(tmp.coef) # coefficients at the right scale including beta0
	structure( c(log(tmp.coef[1]), tmp.coef[2], log(tmp.coef[3])), names=c("beta0l","beta1","beta2l") )
}



#coef.kinRespRepBetaPos <- function(tmp.m){ #coef is used by confint and others
R.methodsS3::setMethodS3("coefKinrespBeta","default", function( 
	### Transform coefficients from Beta-Model from log-Scale to original scale.
	tmp.coef	## beta coefficients at transformed scale, e.g. \code{coef(model1)}
	,...
){
	##details<< 
	## Coefficients beta0 and beta2 are not fitted directyl, but their log is fitted.
	## This ensures that their back-transformation is log-normally distributed and boundbed to be strictly positive.
	
	##details<<
	## If tmp.coef does not include parameter beta0, it is included with default zero.
	
	##seealso<< 
	## \code{\link{fitKinrespBetaReplicate}}
	if( names(tmp.coef)[1] == "beta1" ){
		tmp.coef <- c( beta0=0, tmp.coef )
	}
	if( names(tmp.coef)[1]=="beta0l"){
		tmp.coef <- structure( c(exp(tmp.coef[1]), tmp.coef[2:3]), names=c("beta0",names(tmp.coef)[2:3]) )
	}
	if( "beta2l" %in% names(tmp.coef)){
		tmp.names <- names(tmp.coef)
		tmp.names[match( "beta2l", names(tmp.coef))] <- "beta2"
		names(tmp.coef) <- tmp.names
		tmp.coef["beta2"] <- exp(tmp.coef["beta2"])
	}
	tmp.coef
	### Named vector beta0,beta1,beta2 at original scale.
})

R.methodsS3::setMethodS3("coefKinrespBeta","kinrespList", function(
		### Check the microbial coefficients and translate to original microbial scale for all replicates.
		tmp.coef	##<< result of \code{\link{kinrespGrowthphaseExperiment}}
		,rds.e=NULL ##<< constrained dataset, which may omit some replicates
		,...
	){
		# coefKinresp.kinrespList
		##seealso<< 
		## \code{\link{coefKinresp.default}}
		## ,\code{\link{twKinresp}}
		
		#resRepI <- tmp.coef$resRep[[1]]
		bo <- if( is.null(rds.e)) TRUE else{
				serUnique <- unique(getSERId(rds.e))
				names(tmp.coef$resRep) %in% serUnique
			}		
		tmp <- lapply( tmp.coef$resRep[bo], function(resRepI){ as.data.frame(c(list( experiment=resRepI$dataset$experiment[1], replicate=resRepI$dataset$replicate[1]), coefKinrespBeta.default(coef(resRepI$fit)) ))})
		do.call("rbind",tmp)
		### named numer matrix (columns experiment, replicate, mumax, x0, and r0) with rows corresponding replicates
	})

R.methodsS3::setMethodS3("coefList","kinrespList", function(
		### Extract the coefficients for each replicate
		tmp.coef	##<< result of \code{\link{kinrespGrowthphaseExperiment}}
		,rds.e=NULL ##<< constrained dataset, which may omit some replicates
		,...
	){
		# coefKinresp.kinrespList
		##seealso<< 
		## \code{\link{coefKinresp.default}}
		## ,\code{\link{twKinresp}}
		
		bo <- if( is.null(rds.e)) TRUE else{
				serUnique <- unique(getSERId(rds.e))
				names(tmp.coef$resRep) %in% serUnique
			}		
		#resRepI <- tmp.coef$resRep[[1]]
		tmp <- lapply( tmp.coef$resRep[bo], function(resRepI){ as.data.frame(c(list( experiment=resRepI$dataset$experiment[1], replicate=resRepI$dataset$replicate[1]), (coef(resRepI$fit)) ))})
		do.call("rbind",tmp)
		#tmp
		### named numer matrix (columns experiment, replicate, mumax, x0, and r0) with rows corresponding replicates
	})



.confintKinrespBeta <- function(
	### Transform confidence interval from log-Scale to original scale.
	tmp.cf
){
	if( rownames(tmp.cf)[1]=="beta0l"){
		tmp.cf[1,] <- exp(tmp.cf[1,])
	}
	if( "beta2l" %in% rownames(tmp.cf)){
		tmp.rownames <- rownames(tmp.cf)
		tmp.rownames[match( "beta2l", rownames(tmp.cf))] <- "beta2"
		rownames(tmp.cf) <- tmp.rownames
		tmp.cf["beta2",] <- exp(tmp.cf["beta2",])
	}
	if( rownames(tmp.cf)[1] == "beta1" ){
		tmp.cf <- rbind( c( 0,0), tmp.cf )
	}
	rownames(tmp.cf) <- paste("beta",0:2,sep="")
	tmp.cf
}

.tmp.f <- function(rder.e){
	x <- rder.e$time
	y <- rder.e$resp
}



#-------------- calculating growth parameters from beta fit -----------
calcKinrespCoef <- function(
	### Calculating microbial parameters from beta fit.
	tmp.coef		##<< named numeric vector "beta0" \dots "beta2"
	, lambda=0.9	##<< Ratio of growth associated (coupled) specific respiration to total specific respiration. Usually 0.9.
	, YCO2=1.5		##<< Ratio of assimilated carbon per respired carbon.  Usually 1.5.
	, cf95=NA		##<< Confidence intervals of the transformed beta coefficients. If they are given as numeric matrix with one name row per coefficient, then also cf of the growth parameters are calculated.
){
	##seealso<< 
	## \code{\link{fitKinrespBetaReplicate}}
	
	# if the fit was done with beta0 contrained to 0, we need to add it to tmp.coef
	# and account for log-fit
	tmp.coef <- coefKinrespBeta(tmp.coef)
	
	mumax <- tmp.coef["beta2"]
	r0 <- tmp.coef["beta1"]*(1-lambda) / (tmp.coef["beta0"] +tmp.coef["beta1"]*(1-lambda))
	Q <- mumax / (lambda*YCO2)
	x0 <- tmp.coef["beta1"]/(r0*Q) 
	#x0t <- YCO2*(tmp.coef["beta0"]+tmp.coef["beta1"]*(1-lambda))/((1/lambda-1)*tmp.coef["beta2"])
	if( (x0 < 0) | (r0>1) ){
		res <- structure(c(NA, NA, NA), names=c("mumax","r0","x0"))
	}else{
		res <- structure(c(mumax, r0, x0), names=c("mumax","r0","x0"))
	}
	if( !any(is.na(cf95)) ){
		#convert to original scale involving x0
		tmp <- cf95 <- .confintKinrespBeta(cf95)
		rownames(tmp) <- c("mumax","r0","x0")
		tmp["mumax",] <- cf95["beta2",]
		tmp["r0",1] <- cf95["beta1",1]*(1-lambda) / (cf95["beta0",2] +cf95["beta1",1]*(1-lambda))
		tmp["r0",2] <- cf95["beta1",2]*(1-lambda) / (cf95["beta0",1] +cf95["beta1",2]*(1-lambda))
		tmp["x0",1] <- YCO2*(cf95["beta0",1]+cf95["beta1",1]*(1-lambda))/((1/lambda-1)*cf95["beta2",2]) 
		tmp["x0",2] <- YCO2*(cf95["beta0",2]+cf95["beta1",2]*(1-lambda))/((1/lambda-1)*cf95["beta2",1]) 
		attr(res,"cf95") <- tmp
	}
	res
}

fitKinrespBetaSuite <- function(
	### Fitting the beta-model to all experiments with fixed effects.
	rds.e			##<< dataset with columns experiment, replicate, resp and time, containing only unlimited growth phase (see \code{\link{getUnlimitedGrowthData.kinrespList}}).
	, repFits		##<< Initial coefficients beta0l,beta1, and beta2l for each replicate 
	, weights=NULL	##<< Variance function. see details
){
	if( length(unique(rds.e$suite)) != 1)
		stop("fitKinrespBetaSuite: found other than 1 unique suite identifier in argument rds.e")
	
	rds.eg <- groupedData( resp ~ time | experiment/replicate, data=cbind(rds.e, exprep=(rds.e$experiment:rds.e$replicate)[drop=TRUE] ))
	iFinite <- which(!apply( repFits[,c("beta0","beta2")], 1, function(row){ any(row <= 0) }))
	if( length(iFinite) < nrow(repFits) ){
		repFits <- repFits[iFinite,]
		ser <- getSERId(rds.eg)
		rds.eg <- rds.eg[ser %in% names(iFinite),]
	}
	tmp.exp <- (repFits$experiment)[drop=TRUE] 		
	startRep <- c( log(repFits[,"beta0"]), repFits[,"beta1"],log(repFits[,"beta2"]) )
	# fixed effects for all coefficients by experiment
	tmp.fit <- try(gnls( resp ~ exp(beta0l) + beta1 * exp( exp(beta2l) * time)
			, params=list(beta0l+beta1+beta2l~exprep)
			, start=startRep
			, weights=weights
			, data=rds.eg
		))
	tmp.fit
}

