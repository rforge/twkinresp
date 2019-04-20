setMethodS3("kinrespParDist","gnls", function(
	### Expected value and confidence interval of microbial parameters
	model	##<< the result of the \code{\link{fitKinrespReplicate}} or \code{\link{fitKinrespExperiment}$model}
	,...	##<< currently not used
){
	# kinrespParDist.gnls
	##alias<< kinrespParDist
	##seealso<<
	## \code{\link{coefKinresp.default}}
	## ,\code{\link{twKinresp}}

	##details<<
	## Both the lognormal and the logitnormal distributions are skewed.
	## Therefore, the mode, the median, and the expected value differ.
	## At transformed, i.e. normal scale, they coincide. Simple backtransformation
	## yield the median at original scale.

	#mu0 <- coef(model)
	mu0 <- fixef(model)		#works on both gnls and nlme
	if (!all(c("mumaxl","r0l","x0l") %in% names(mu0)) )
		stop("kinrespParDist.gnls: provided model with some normal scale estimates missing.")
	sigma0 <- summary(model)$tTable[,"Std.Error"]	#same as sqrt(diag(model$varBeta)) for gls
	sigma02 <- sigma0^2
	cf0 <- confint(model)
	# first assume lognormal distribution
	colNames <- c("mean","sd","mle","median","cf025","cf975","mu","sigma")
	res <- matrix(as.numeric(NA), nrow=3, ncol=length(colNames), dimnames=list( c("x0","r0","mumax"), colNames ))
	#log-normal
	i0 <- match(c("x0l","mumaxl"), names(mu0) )
	mu <- mu0[i0]; sigma <- sigma0[i0]; sigma2 <- sigma02[i0]; cf<-cf0[i0,]
	i <- match(c("x0","mumax"), rownames(res) )
	res[i,"mean"] <- exp(mu+sigma2/2)
	res[i,"sd"] <- sqrt( (exp(sigma2)-1)*exp( 2*mu+sigma2 ) )
	res[i,"mle"] <- exp(mu-sigma2)
	res[i,"median"] <- exp(mu)
	res[i,c("cf025","cf975")] <- exp(cf)
	res[i,c("mu","sigma")] <- cbind(mu,sigma)
	#logit-normal
	i0 <- match(c("r0l"), names(mu0) )
	mu <- mu0[i0]; sigma <- sigma0[i0]; sigma2 <- sigma02[i0]; cf<-cf0[i0,]
	i <- match(c("r0"), rownames(res) )
	moments <- momentsLogitnorm( mu, sigma)
	res[i,"mean"] <- moments["mean"]
	res[i,"sd"] <- sqrt(moments["var"])
	#mtrace(modeLogitnorm)
	res[i,"mle"] <- modeLogitnorm(mu,sigma)
	res[i,"median"] <- invlogit(mu)
	res[i,c("cf025","cf975")] <- invlogit(cf)
	res[i,c("mu","sigma")] <- cbind(mu,sigma)
	res
	### numeric matrix with rows corresponding to variables and columns \itemize{
	### \item{ \code{mean}: expected value}
	### \item{ \code{sd}: standard deviation}
	### \item{ \code{mle}: the mode, i.e. the maximum likelihood estimate}
	### \item{ \code{median}: the median}
	### \item{ \code{cf025} and \code{cf975}: 2.5% and 97.5% confidence bounds}
	### \item{ \code{mu} and \code{sigma}: parameters at normal, i.e log or logit transformed scale}
	### }
})

setMethodS3("kinrespParDist","nlme", function(
	### Expected value and confidence interval of microbial parameters
	model	##<< the result of the \code{\link{fitKinrespReplicate}} or \code{\link{fitKinrespExperiment}$model}
	,...	##<< currently not used
){
	##seealso<<
	## \code{\link{kinrespParDist.gnls}}
	## ,\code{\link{coefKinresp.default}}
	## ,\code{\link{twKinresp}}
	kinrespParDist.gnls(model,...)
})

setMethodS3("coefKinresp","default", function(
		### Check the microbial coefficients and translate to original microbial scale.
		tmp.coef	##<< coefficients of model fit \code{coef(model)} (see details)
		,...
	){
		# coefKinresp.default
		##alias<< coefKinresp
		##seealso<<
		## \code{\link{twKinresp}}

		##details<< \describe{\item{Functions for acccessing microbial parameters, and their uncertainty bounds.}{
		## \itemize{
		## \item{ Mode, Median, Mean and confidence bounds of microbial fits: \code{\link{kinrespParDist.gnls}}}
		## \item{ Microbial parameters of single fit: this method (obtained by \code{coef(\link{kinrespGrowthphaseReplicate})} or \code{fixef(\link{fitKinrespExperiment})}) }
		## \item{ Microbial parameters of replicate fits: \code{\link{coefKinrespMatrix}} (obtained by \code{coef(\link{fitKinrespExperiment})}) }
		## \item{ Microbial parameters of a list of fits: \code{\link{coefKinresp.kinrespList}} (obtained by \code{\link{kinrespGrowthphaseExperiment}})}
		## \item{ 95% confidenc interval of microbial parameters of a single fit: \code{\link{confintKinresp}} }
		## \item{ Microbial parameters from beta-form of model fit: \code{\link{calcKinrespCoef}} }
		## }
		##}}

		##details<< \describe{\item{Functions for translating between normal and original scale.}{
		## \itemize{
		## \item{ normalized to original scale: all the methods above. }
		## \item{ original scale to normalized scale: \code{\link{coefKinrespNormStart}} }
		## }
		##}}

		##details<< \describe{\item{ \code{coef(model)} }{
		## Models are fitted in various forms differing by used coefficients.
		## This method recognizes and translates coefficients of the following forms
		## \itemize{
		## \item{ beta-form at original scale and at normalized scale}
		## \item{ beta-form excluding beta0 (fixed to 0 see
		## \code{\link{calcKinrespCoef}}) }
		## \item{ microbial form fit at original and transformed scale.}
		## \item{ microbial form fit with excluding r0 (assuming r0=1)}
		## }
		## }}
		if (names(tmp.coef)[1] %in% c("beta0l","beta0","beta1")) {
			tmp.coef <- calcKinrespCoef(tmp.coef)
		}
		if ("x0l" %in% names(tmp.coef)) {
			tmp.names <- names(tmp.coef)
			tmp.names[match( "x0l", names(tmp.coef))] <- "x0"
			names(tmp.coef) <- tmp.names
			tmp.coef["x0"] <- exp(tmp.coef["x0"])
		}
		if ("r0l" %in% names(tmp.coef)){
			tmp.names <- names(tmp.coef)
			tmp.names[match( "r0l", names(tmp.coef))] <- "r0"
			names(tmp.coef) <- tmp.names
			tmp.coef["r0"] <- invlogit(tmp.coef["r0"])
		}
		if ("mumaxl" %in% names(tmp.coef)){
			tmp.names <- names(tmp.coef)
			tmp.names[match( "mumaxl", names(tmp.coef))] <- "mumax"
			names(tmp.coef) <- tmp.names
			tmp.coef["mumax"] <- exp(tmp.coef["mumax"])
		}
		if (!("r0" %in% names(tmp.coef)) ){
			tmp.r0 = structure( attr(tmp.coef,"r0"), names=NULL )
			if (is.null(tmp.r0) ) tmp.r0 = 1
			tmp.coef <- c(tmp.coef["mumax"],r0=tmp.r0,tmp.coef["x0"])
		}
		tmp.coef
		### named numeric vector (x0,r0,mumax)
	})

setMethodS3("coefKinresp","kinrespList", function(
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
		bo <- if (is.null(rds.e)) TRUE else{
				serUnique <- unique(getSERId(rds.e))
				names(tmp.coef$resRep) %in% serUnique
			}
		tmp <- lapply( tmp.coef$resRep[bo], function(resRepI){ as.data.frame(c(list( experiment=resRepI$dataset$experiment[1], replicate=resRepI$dataset$replicate[1]), coefKinresp.default(coef(resRepI$fit)) ))})
		do.call("rbind",tmp)
		### named numer matrix (columns experiment, replicate, mumax, x0, and r0) with rows corresponding replicates
	})

setMethodS3("fixef","kinresp", function(
		### Make fixed.effects deliver attribute r0 if this is supplied with model.
		object
		,...
	){
		##seealso<<
		## \code{\link{coefKinresp.default}}
		## ,\code{\link{twKinresp}}
		attr(object,"class") <- class(object)[-1]
		tmp <- fixef(object)
		if (!is.null(attr(object,"r0")) )
			attr(tmp,"r0") <- attr(object,"r0")
		tmp
	})

setMethodS3("coef","kinresp", function(
		### Make coef deliver attribute r0 if this is supplied with model.
		object
		,...
	){
		# coef.kinresp
		##seealso<<
		## \code{\link{coefKinresp.default}}
		## ,\code{\link{twKinresp}}
		attr(object,"class") <- class(object)[-1]
		tmp <- coef(object)
		if (!is.null(attr(object,"r0")) )
			attr(tmp,"r0") <- attr(object,"r0")
		tmp
	})

setMethodS3("confint","kinresp", function(
		### Make confint deliver attribute r0 if this was supplied with object
		object
		,...
	){
		##seealso<<
		## \code{\link{coefKinresp.default}}
		## ,\code{\link{twKinresp}}
		attr(object,"class") <- class(object)[-1]
		tmp <- confint(object)
		if (!is.null(attr(object,"r0")) )
			attr(tmp,"r0") <- attr(object,"r0")
		tmp
	})

coefKinrespMatrix <- function(
	### Check the microbial coefficients and translates to original microbial scale.
	tmp.coef	##<< multiple rows of microbial coefficients (see \code{\link{coefKinresp.default}})
){
	##seealso<<
	## \code{\link{coefKinresp.default}}
	## ,\code{\link{twKinresp}}
	if (colnames(tmp.coef)[1] %in% c("beta0l","beta0","beta1") ){
		tmp.coef <- calcKinrespCoef(tmp.coef)
	}
	if ("x0l" %in% colnames(tmp.coef)){
		tmp.colnames <- colnames(tmp.coef)
		tmp.colnames[match( "x0l", colnames(tmp.coef))] <- "x0"
		colnames(tmp.coef) <- tmp.colnames
		tmp.coef[,"x0"] <- exp(tmp.coef[,"x0"])
	}
	if ("r0l" %in% colnames(tmp.coef)){
		tmp.colnames <- colnames(tmp.coef)
		tmp.colnames[match( "r0l", colnames(tmp.coef))] <- "r0"
		colnames(tmp.coef) <- tmp.colnames
		tmp.coef[,"r0"] <- invlogit(tmp.coef[,"r0"])
	}
	if ("mumaxl" %in% colnames(tmp.coef)){
		tmp.colnames <- colnames(tmp.coef)
		tmp.colnames[match( "mumaxl", colnames(tmp.coef))] <- "mumax"
		colnames(tmp.coef) <- tmp.colnames
		tmp.coef[,"mumax"] <- exp(tmp.coef[,"mumax"])
	}
	if (!("r0" %in% colnames(tmp.coef)) ){
		tmp.r0 = structure( attr(tmp.coef,"r0"), colnames=NULL )
		if (is.null(tmp.r0) ) tmp.r0 = 1
		tmp.coef <- c(tmp.coef["mumax"],r0=tmp.r0,tmp.coef["x0"])
	}
	tmp.coef
}

confintKinresp <- function(
	### Translate the confidence interval of estimated coefficients to original scale.
	tmp.cf		##<< confidence interval from fitting the microbial model \code{confint(modelfit)}

){
	# confintKinresp
	##seealso<<
	## \code{\link{coefKinresp.default}}
	## ,\code{\link{twKinresp}}
	tmp.i <- match("x0l", rownames(tmp.cf))
	if (!is.na(tmp.i)){
		tmp.cf[tmp.i,] <- exp(tmp.cf[tmp.i,])
		rownames(tmp.cf)[tmp.i] <- "x0"
	}
	tmp.i <- match("r0l", rownames(tmp.cf))
	if (!is.na(tmp.i)){
		tmp.cf[tmp.i,] <- invlogit(tmp.cf[tmp.i,])
		rownames(tmp.cf)[tmp.i] <- "r0"
	}
	tmp.i <- match("mumaxl", rownames(tmp.cf))
	if (!is.na(tmp.i)){
		tmp.cf[tmp.i,] <- exp(tmp.cf[tmp.i,])
		rownames(tmp.cf)[tmp.i] <- "mumax"
	}
	if (rownames(tmp.cf)[2] != "r0" ){
		tmp.r0 = structure( attr(tmp.cf,"r0"), names=NULL )
		if (is.null(tmp.r0) ) tmp.r0 = 1
		tmp.cf <- rbind( tmp.cf["mumax",],c( tmp.r0,tmp.r0), tmp.cf["x0",] )
		rownames(tmp.cf) <- c("mumax","r0","x0")
	}
	tmp.cf
}

.coefKinrespLogStart <- function(
	### Transform microbial coefficients to log scale (does ot care for r0)
	tmp.coef
){
	##seealso<<
	## \code{\link{coefKinresp.default}}
	## ,\code{\link{twKinresp}}

	##seealso<< \code{\link{coefKinrespNormStart}}
	tmp.coef <- coefKinresp(tmp.coef) # coefficients at the original scale including x0
	tmp.names <- names(tmp.coef)
	tmp.names[match( "x0", names(tmp.coef))] <- "x0l"
	tmp.names[match( "mumax", names(tmp.coef))] <- "mumaxl"
	names(tmp.coef) <- tmp.names
	tmp.coef["x0l"] <- log(tmp.coef["x0l"])
	tmp.coef["mumaxl"] <- log(tmp.coef["mumaxl"])
	tmp.coef
}

coefKinrespNormStart <- function(
	### Transform microbial coefficients to normal scale.
	tmp.coef	##<< starting parameters of kinetic respiration, supplied to \code{\link{coefKinresp.default}} (maybe beta)
){
	##seealso<<
	## \code{\link{coefKinresp.default}}
	## ,\code{\link{twKinresp}}

	# coefKinrespNormStart
	##details<<
	## Replaces \itemize{
	## \item{x0 by x0l = log(x0)}
	## \item{mumax by mumaxl=log(mumax)} and
	## \item{r0 by r0l=logit(r0l)}
	## }
	## Values at normalized scale are used as starting values for model fits.
	tmp.coef <- coefKinresp(tmp.coef)
	tmp.names <- names(tmp.coef)
	tmp.names[match( "x0", names(tmp.coef))] <- "x0l"
	tmp.names[match( "r0", names(tmp.coef))] <- "r0l"
	tmp.names[match( "mumax", names(tmp.coef))] <- "mumaxl"
	names(tmp.coef) <- tmp.names
	tmp.coef["x0l"] <- log(tmp.coef["x0l"])
	tmp.coef["r0l"] <- logit(tmp.coef["r0l"])
	tmp.coef["mumaxl"] <- log(tmp.coef["mumaxl"])
	tmp.coef
	### named vector of coefficients (x0l,r0l,mumaxl)

}

