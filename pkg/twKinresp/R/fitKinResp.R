# function for fitting a general three parameter exponential model
# and functions for fitting the kinRespModel with microbial parameters directly

#----------------- fitting microbial parameters directly -----------------
# see fitKinresp_test.R
fitKinrespExperiment <- function(
	### Fit microbial form to several replicated respiration time series of one experiment.
	rde.e			##<< dataset with columns experiment, replicate, resp and time, containing only unlimited growth phase (see \code{\link{getUnlimitedGrowthData.kinrespList}}).
	, repFits		##<< Initial coefficients mumax, x0, and r0 for each replicate (see \code{\link{coefKinresp.kinrespList}}).
	, lambda=0.9	##<< Ratio of growth associated (coupled) specific respiration to total specific respiration. Usually 0.9.
	, YCO2=1.5		##<< Ratio of assimilated carbon per respired carbon (Y/(1-Y)).  Usually 1.5.
	, weights=NULL	##<< Variance function. see details
	, start=NULL	##<< Specifying starting values in an alternative way to argument repFits. 
	, tmp.names = c("none","x0l","r0l","x0l+r0l","x0l+r0l+mumaxl")	##<< scenarios of random effects
	, showFitErrorMsg=FALSE	##<< if FALSE (standard) then error msg are suppressed when no fit is obtained. This is a common case for the test of variants that include too many random effects.
){
    ##details<< \describe{\item{Microbial parameters across replicates}{
    ## Microbial parameters can be inferred from fitting against the time series of respiration 
    ## of each replicate separately (\code{\link{fitKinrespReplicate}}.
    ## However, what are then parameters of the population? The average across replicate parameters will be wrong.
    ##
    ## Alternatively, one can first average the measurement across replicates at for each measurement time.
    ## However, the uncertainty of the parameters will be wrong.
    ##
    ## A viable soluion is to fit a mixed model to all the replicate data. The micoribal parameters, then are
    ## described by a mean value of the population and a variance across replicates.
    ## Several options of which parameters vary across replicates can be tested.
    ## }}
    
    ##details<<
    ## \describe{\item{Variance function}{
	## If no weights are given, the measurement errors of the single
	## observations are assumed to be identical. If measurement errors
	## increase with the magnitude of the observations, this can be
	## modeled by power variance function:
	## \code{weights=\link{varPower}(fixed=delta)},
	## with delta being a value between 0 (constant expected absolute error) 
	## and 1 (constant expected relative error).
	## }}
	
	if( length(unique(rde.e$experiment)) != 1)
		stop("fitKinrespExperiment: found other than 1 unique experiment identifier in argument rde.e")
	
	if( is.null(start )){
		if( is.null(repFits) ) stop("fitKinExperimentSuite: must provide starting values (start) or estimates for each replicate (repFits).")
		tmp.coefs <- subset(repFits, repFits$experiment %in% unique(rde.e$experiment))
		tmp.mumaxl <- mean(log(tmp.coefs[,"mumax"]), na.rm=TRUE)
		tmp.x0l <- mean(log(tmp.coefs[,"x0"]), na.rm=TRUE)
		tmp.r0l <- mean(logit(tmp.coefs[,"r0"]), na.rm=TRUE)
		start = c( tmp.mumaxl, tmp.r0l, tmp.x0l  ) # r0 comes before x0 in coef
	}
	
	#the same rep might occur in several experiment, but random factor should treat it as differnt, create unqiue rep
	#better use grouped factors 
	#rde.e$exprep = paste( rde.e$experiment, rde.e$replicate, sep="_")
	rde.e$exprep <- with(rde.e, experiment:replicate)[drop=TRUE]
	rde.eg <- groupedData( resp ~ time | experiment/replicate, data=cbind(rde.e, lambda=lambda, YCO2=YCO2))
	# as a parameter tmp.names <- c("none","x0l","r0","x0l+r0","x0l+r0+mumaxl")
	tmp.AIC <- structure( rep(NA, length(tmp.names)),names=tmp.names)
	tmp.fits <- list()
	try(suppressWarnings({
				tmp.name <- tmp.names[1]
				tmp.fit1 <- gnls(
					resp ~ exp(x0l)*(1-invlogit(r0l))*(1/lambda-1)*exp(mumaxl)/YCO2 + exp(x0l)*invlogit(r0l)*1/lambda*exp(mumaxl)/YCO2 * exp( exp(mumaxl) * time)  
					, params=list(mumaxl+r0l+x0l~1)
					, start=start
					, weights=weights
					, data=rde.eg
				)
				tmp.fits[[tmp.name]] <- tmp.fit1
				tmp.AIC[tmp.name] <- AIC(tmp.fit1)
			}),silent=!showFitErrorMsg)
	#try various models of random effects in different coefficients
	tmp.name <- "x0l"
	for( tmp.name in tmp.names[-1]){
		#tmp.random <- eval(parse(text=paste(tmp.name,"~1|experiment/replicate",sep="")))
		tmp.random <- eval(parse(text=paste(tmp.name,"~1|exprep",sep="")))
		# not accepted as random formula: tmp.random <- eval(parse(text=paste(tmp.name,"~1|experiment:replicate",sep="")))
		#tmp.random <- eval(parse(text=paste("list(experiment=NULL, exprep=",tmp.name,"~1)",sep="")))
		#tmp.random <- eval(parse(text=paste("list(experiment=",1,"~1, exprep=",tmp.name,"~1)",sep="")))
		#tmp.random <- eval(parse(text=paste("list(experiment=NULL, exprep=",tmp.name,"~1)",sep="")))
		if( showFitErrorMsg ) print(tmp.random)
		try(suppressWarnings({
					tmp.fit2 <- nlme(
						resp ~ exp(x0l)*(1-invlogit(r0l))*(1/lambda-1)*exp(mumaxl)/YCO2 + exp(x0l)*invlogit(r0l)*1/lambda*exp(mumaxl)/YCO2 * exp( exp(mumaxl) * time)  
						#,fixed=list(mumaxl+r0l+x0l~experiment)	#include fixed effects for each experiment
						,fixed=list(mumaxl+r0l+x0l~1)	#include fixed effects for each experiment
						,random=tmp.random
						, weights =weights
						, start=start
						, method='REML'	# for unbiased estimates of standard error and confidence intervals 
						,data=rde.eg
					)
					tmp.fits[[tmp.name]] <- tmp.fit2
					tmp.AIC[tmp.name] <- AIC(tmp.fit2)
				}),silent=!showFitErrorMsg)
	}
	tmp.random <- tmp.names[ which.min(tmp.AIC)]
	tmp.fit <- tmp.fits[[ tmp.random  ]]
	list( 
		model=tmp.fit			
		, random=tmp.random		
		, fits=tmp.fits			
		, aics=tmp.AIC 			 
	)
	### list with components \describe{
	### \item{model}{best fitting result}
	### \item{random}{best random effects scenario}
	### \item{fits}{fits of all random effects scenarios}
	### \item{aics}{Akaike information criterion for random effects scenarios}
	### }
}
attr(fitKinrespExperiment,"ex") <- function(){
    # data of one example treament: measurements of several replicates of one soil
	data(respWutzler10)
	rde <- subset(respWutzler10, suite=="Face" & experiment==9 )
	
	# constrain data to unlimited growth phase
	res4 <- kinrespGrowthphaseExperiment(rde, weights=varPower(fixed=0.5) )
    rde.e <- getUnlimitedGrowthData(res4)
	# fit the mixed model to all replicates
	res5Scen <- fitKinrespExperiment( rde.e, coefKinresp(res4,rde.e), weights=varPower(fixed=0.5) )
	
    # get the best fit parameters of the population
    coefKinresp(fixef(res5Scen$model))
    
	# examine the random-effects scenarios: 
	# Here the lowest AIC suggest that activity state and initial microbial biomass
	# differed between replicates, but maximum growth was the same
	res5Scen$aics
	
	# plot the fits
	#windows(record=TRUE)
	rde.e$fitted <- fitted(res5Scen$model)
	plot( resp ~ time, data=rde.e, col=rde.e$replicate )
	tmp <- by( rde.e, rde.e$replicate, function(rder){ lines(fitted~time,data=rder, col=as.numeric(rder[1,replicate]))})
	
	# estimated microbial coefficients
	(pars <- kinrespParDist(res5Scen$model))
	# plot the density of r0 and density summaries
	iPar="r0"
	xGrid <- seq( pars[iPar,"cf025"]*0.8, pars[iPar,"cf975"]*1.2, length.out=80)
	#fx <- dlnorm(xGrid, mean=pars[iPar,"mu"],sd=pars[iPar,"sigma"])
	fx <- dlogitnorm(xGrid, mu=pars[iPar,"mu"],sigma=pars[iPar,"sigma"])
	plot( fx ~ xGrid, type="l", xlab=iPar, ylab="density" )
	abline(v=pars[iPar,c("mle","median","mean","cf025","cf975")], col=c("red","green","blue","gray","gray"))
}

modelKinrespMic <- function(
	### Calculate respiration for time, using microbial parameters 
	time	
	,param			##<< named numeric vector ("x0","r0","mumax")
	, lambda=0.9	##<< Ratio of growth associated (coupled) specific respiration to total specific respiration. Usually 0.9.
	, YCO2=1.5		##<< Ratio of assimilated carbon per respired carbon.  Usually 1.5.
){
	##seealso<< 
	## \code{\link{twKinresp}}
	
	param["x0"]*(1-param["r0"])*(1/lambda-1)*param["mumax"]/YCO2 + param["x0"]*param["r0"]*1/lambda*param["mumax"]/YCO2 * exp( param["mumax"] * time)
	### \code{x0*(1-r0)*(1/lambda-1)*mumax/YCO2 + x0*r0*1/lambda*mumax/YCO2 * exp(mumax*time)}
}

fitKinrespReplicate <- function(
	### Fitting the kinResp microbial model to time series of one replicate.
	rder.e			##<< Respiration dataset containing columns replicate, resp and time, constrained to unlimited growth phase.
	, lambda=0.9	
	, YCO2=1.5
	, start=NULL
	, weights=NULL 
){
	# fitKinrespRepl

	##seealso<< 
	## \code{\link{twKinresp}}, \code{\link{coefKinresp}} 
	
	##details<<  
	## If the microbial explicit form did not fit, then the beta fit is returned.
	## The beta fit itself might not contain beta0, because it was constrained to 1.
	
	if( length(unique(rder.e$replicate)) != 1)
		stop("fitKinrespExperiment: found other than 1 unique replicate identifier in argument rder.e")
	
	if( is.null(start) ){
		#obtain the starting values from linear fit of the beta model form
		start.beta <- coef(lm(log(rder.e$resp-0.99*min(rder.e$resp)) ~ rder.e$time))
		tmp.beta=structure(
			c(beta0=0.99*min(rder.e$resp),beta1=exp(start.beta[1]),beta2=start.beta[2])
			,names=paste("beta",0:2,sep="")
		)
		#tmp.beta <- coef(fitExpModel(x=rder.e$time, y=rder.e$resp))
		start <- calcKinrespCoef(tmp.beta, lambda=lambda, YCO2=YCO2)
	}
	# replace starting values by logarithm
	#start.log <- coefKinrespLogStart(start)
	start.norm <- coefKinrespNormStart(start)
	tmp.fit <- try({
			#tmp.fit <- gnls( resp ~ exp(x0l)*(1-r0)*(1/lambda-1)*mumax/YCO2 + exp(x0l)*r0*1/lambda*mumax/YCO2 * exp( mumax * time)
			tmp.fit <- gnls( resp ~ exp(x0l)*(1-invlogit(r0l))*(1/lambda-1)*exp(mumaxl)/YCO2 + exp(x0l)*invlogit(r0l)*1/lambda*exp(mumaxl)/YCO2 * exp( exp(mumaxl) * time)
				#, params=mumax+r0+x0l~1
				, params=mumaxl+r0l+x0l~1
				, start=start.norm
				, weights=weights
				, data=cbind(rder.e, lambda=lambda, YCO2=YCO2)
			)
			#tmp.r0 <- logit(coef(tmp.fit)["r0l"])
			tmp.fit
		})
	# fall back one, use start parameters of beta-Fit
	if( inherits(tmp.fit,"try-error" ) ){
		tmp.beta <- coefKinrespBeta(coef(tmp <- fitKinrespBetaReplicate(x=rder.e$time, y=rder.e$resp)))
		start <- calcKinrespCoef(tmp.beta, lambda=lambda, YCO2=YCO2, cf95=confint(tmp))
		#start.log <- coefKinrespLogStart(start)
		start.norm <- coefKinrespNormStart(start)
		tmp.fit <- try(
			if( start.norm["r0"] == 1){
					# if r0 > 1 refit with constraining r0 to 1
					tmp.fit <- gnls( resp ~ exp(x0l)*1/lambda*exp(mumaxl)/YCO2 * exp( exp(mumaxl) * time)
						, params=mumaxl+x0l~1	#without r0 that will be fixed to 1
						, start=start.norm[c(1,3)]
						, weights=weights
						, data=cbind(rder.e, lambda=lambda, YCO2=YCO2, r0=1)
					)
					attr(tmp.fit,"r0") <- 1
					attr(tmp.fit,"class") <- c("kinresp",class(tmp.fit))
					tmp.fit
				}else{ 
					gnls( resp ~ exp(x0l)*(1-invlogit(r0l))*(1/lambda-1)*exp(mumaxl)/YCO2 + exp(x0l)*invlogit(r0l)*1/lambda*exp(mumaxl)/YCO2 * exp( exp(mumaxl) * time)
						, params=mumaxl+r0l+x0l~1
						, start=start.norm
						, weights=weights
						, data=cbind(rder.e, lambda=lambda, YCO2=YCO2)
					)				
				}
		)
		# fall back two, deliver the beta model
		if( inherits(tmp.fit,"try-error" ) ){
			tmp.fit <- tmp
		}
	}
	tmp.fit
}


.tmp.f <- function(rder.e){
	(tmp.m <- fitKinrespReplicate(rder.e))
	plot(rder.e$resp ~ rder.e$time)
	(tmp.mb <- fitExpModel(x=rder.e$time, y=rder.e$resp))
	lines( modelKinrespBeta(rder.e$time,tmp.beta) ~ rder.e$time )
	lines( modelKinrespMic(rder.e$time,tmp.start) ~ rder.e$time, col="red" )
}

#------------ fitting all replicates
.fitKinrespReplAll <- function(
	### Invoking \code{\link{fitKinrespRepl}} for each experiment within the suite.
	rd.e			##<< constrained data series
	, lambda=0.9
	, YCO2=1.5
	, weights=NULL 
){
	# fitKinrespReplAll
	##seealso<< 
	## \code{\link{twKinresp}}
	
	coefRep <- data.frame()
	mRep <- list()
	i_exp <- "Ek_15_27"
	i_exp <- "17"
	for( i_exp in as.character(unique(rd.e$experiment))){
		print(i_exp)
		rde.e <- subset( rd.e, rd.e$experiment==i_exp )
		i_rep <- "22"
		i_rep <- "1"
		for( i_rep in as.character(unique(rde.e$replicate))){
			print(paste("replicate=",i_rep))
			rder.e <- subset( rde.e, replicate==i_rep)
			tmp.m <- fitKinrespReplicate(rder.e, lambda=lambda, YCO2=YCO2, weights=weights)
			mRep[[ paste(i_exp,i_rep,sep="_") ]] <- tmp.m
			tmp.coef <- coefKinresp(coef(tmp.m))
			coefRep <- rbind( data.frame(exp=i_exp, rep=i_rep, as.list(tmp.coef)), coefRep )
		}
	}
	list( m = mRep, coef=coefRep)
	###	list( m identifier of "idExp_idRep", coef: data.frame with each row one replicate,
	### columns exp, rep, \code{\link{coefKinresp(tmp.m)}}
}
#mtrace(fitKinrespReplAll)

#------------- entire experiment: mixed model to beta model ---------
.fitKinExperimentBeta <- function(
	### Fit an exponential equation to the data of the experiment (each series constrained to unlimited growth phase).
	rde.e
	, weights=NULL
){
	##seealso<< 
	## \code{\link{twKinresp}}
	
	tmp.fit.e <- .fitExpModelConstr(rde.e$time,rde.e$resp, weights=weights)
	rde.e <- groupedData( resp ~ time | replicate, data=rde.e)
	start = list(fixed=coef(tmp.fit.e))
	if( is.null(start$fixed$beta0) ){
		#fit beta0 fixed to 0 
		tmp.names <- c("none","beta1","beta2"
				,"beta1+beta2"
		) 
		tmp.aic <- structure( c(AIC(tmp.fit.e),rep(Inf,length(tmp.names)-1)), names= tmp.names)
		tmp.fits <- list("none"=tmp.fit.e)
		tmp.name <- tmp.names[2]
		for( tmp.name in tmp.names[-1] ){
			tmp.random <- eval(parse(text=paste(tmp.name,"~1",sep="")))
			print(tmp.random)
			try({ 
				tmp.fit2.e <-nlme(
					resp ~ beta1 * exp( beta2 * time) 
					,fixed=beta1+beta2~1
					,random=tmp.random
					, start=start
					, method='REML'	# for unbiased estimates of standard error and confidence intervals 
					,data=rde.e
				)
			},silent=TRUE)
			tmp.fits[[tmp.name]] <- tmp.fit2.e
			tmp.aic[tmp.name] <- AIC(tmp.fit2.e)
		}
	}else{ # three parameter model
		startlog <- start; startlog$fixed$beta0 <- NULL; startlog$fixed$beta0l <- log(start$fixed$beta0)
		# for this case allow full correlations of the parameters, do not work with pdDiag etc
		tmp.names <- c("none","beta0l","beta1","beta2"
				,"beta0l+beta1" 
				,"beta0l+beta1+beta2")
		tmp.aic <- structure( c(AIC(tmp.fit.e),rep(Inf,length(tmp.names)-1)), names= tmp.names)
		tmp.fits <- list("none"=tmp.fit.e)
		tmp.name <- tmp.names[2]
		for( tmp.name in tmp.names[-1] ){
			tmp.random <- eval(parse(text=paste(tmp.name,"~1",sep="")))
			print(tmp.random)
			try({
				tmp.fit2.e <- nlme(
					resp ~ log(beta0l) + beta1 * exp( beta2 * time) 
					,fixed=beta0l+beta1+beta2~1
					,random=tmp.random
					, start=	#starting values from gnls
					, method='REML'	# for unbiased estimates of standard error and confidence intervals 
					,data=rde.e
				)
				tmp.fits[[tmp.name]] <- tmp.fit2.e
				tmp.aic[tmp.name] <- AIC(tmp.fit2.e)
			})
		}
		tmp.f <- function(){
			#tmp.names <- c("none","beta0","beta1","beta2","beta0+beta1","beta0+beta1_cor","beta0+beta1_corf","beta0+beta2","beta0+beta1+beta2","beta0+beta1+beta2_cor","beta0+beta1+beta2_corf") 
			tmp.aic <- structure( c(AIC(tmp.fit.e),AIC(tmp.fit2.e),rep(Inf,length(tmp.names)-2)), names= tmp.names)
			tmp.fits <- list("none"=tmp.fit.e, "beta0"=tmp.fit2.e)
			
			try({ tmp.fits$beta1 <- update( tmp.fit2.e, random=beta1~1); tmp.aic["beta1"] <- AIC(tmp.fits$beta1) },silent = TRUE)	
			try({ tmp.fits$beta2 <- update( tmp.fit2.e, random=beta2~1); tmp.aic["beta2"] <- AIC(tmp.fits$beta2) },silent = TRUE)	
			tmp.fit <- tmp.fits[[names(tmp.aic)[ which.min(tmp.aic)] ]]
			
			#from now in increasing order so that most complex sucessful fit overwrites the previous one
			
			#try({ tmp.fits[["beta0+beta2"]] <- tmp.fit <- update( tmp.fit2.e, random=pdDiag(beta0+beta2~1)); tmp.aic["beta0+beta2"] <- AIC(tmp.fit) },silent = TRUE)
			#try({ tmp.fits[["beta0+beta2_cor"]] <- tmp.fit <- update( tmp.fit2.e, random=pdSymm(beta0+beta2~1)); tmp.aic["beta0+beta2_cor"] <- AIC(tmp.fit) },silent = TRUE)
			#try({ tmp.fits[["beta1+beta2"]] <- tmp.fit <- update( tmp.fit2.e, random=pdDiag(beta1+beta2~1)); tmp.aic["beta1+beta2"] <- AIC(tmp.fit) },silent = TRUE)
			#try({ tmp.fits[["beta1+beta2_cor"]] <- tmp.fit <- update( tmp.fit2.e, random=pdSymm(beta1+beta2~1)); tmp.aic["beta1+beta2_cor"] <- AIC(tmp.fit) },silent = TRUE)
			
			# check on random in initial activity 
			try({ tmp.fits[["beta0+beta1_cor"]] <- tmp.fit <- update( tmp.fit2.e, random=pdSymm(beta0+beta1~1)); tmp.aic["beta0+beta1_cor"] <- AIC(tmp.fit)-2 },silent = TRUE) #correlation parameter of nearly-1 is always due to state
			
			# check on additional effect of  
			try({ 
						#allow correlation between beta0 and beta1 but beta2 independent
						tmp.fits[["beta0+beta1+beta2_cor"]] <- tmp.fit <- update( tmp.fit2.e, random=
										pdBlocked(list(beta0+beta1~1,beta2~1),pdClass=c("pdSymm","pdDiag"))
								) 
						#tmp.fits[["beta0+beta1+beta2_cor"]] <- tmp.fit <- update( tmp.fit2.e, random=pdSymm(beta0+beta1+beta2~1)) 
						tmp.aic["beta0+beta1+beta2_cor"] <- AIC(tmp.fit)
					},silent = TRUE)
		}#tmp.f of more complicated covariance matrices
	}# else omit beta0l
	tmp.random <- names(tmp.aic)[ which.min(tmp.aic)]
	tmp.fit <- tmp.fits[[ tmp.random ]]
	### list with entries
	list( 
		##describe<< 
		fit=tmp.fit	##<< the fit with the lowest AIC may be of class gnls or nlme 
		,fits=tmp.fits	##<< List of results of all the variants of inclusion of random effects. 
		, aics=tmp.aic	##<< Akaike Information criterion for fits,  information on errors, on some replicates
		##end<< 
	)
}
#mtrace(.fitKinExperimentBeta)

.tmp.f <- function(tmp.fit2.e){
	#no sense to fit uncorrelated beta0 and beta1
	try({ tmp.fits[["beta0+beta1"]] <- tmp.fit <- update( tmp.fit2.e, random=pdDiag(beta0+beta1~1)); tmp.aic["beta0+beta1"] <- AIC(tmp.fits[["beta0+beta1"]]) },silent = TRUE)	
	try({ tmp.fits[["beta0+beta1+beta2"]] <- update( tmp.fit2.e, random=pdDiag(beta0+beta1+beta2~1)); tmp.aic["beta0+beta1+beta2"] <- AIC(tmp.fits[["beta0+beta1+beta2"]]) },silent = TRUE)
	#fixing the correlation does not work
	try({
				#allow correlation between beta0 and beta1
				#tmp.fits[["beta0+beta1_cor"]] <- tmp.fit <- update( tmp.fit2.e, random=pdSymm(beta0+beta1~1)); 
				tmp.fits[["beta0+beta1_corf"]] <- tmp.fit <- update( tmp.fits[["beta0+beta1_cor"]], random=pdDiagCorrNeg(beta0+beta1~1)); 
				tmp.aic["beta0+beta1_corf"] <- AIC(tmp.fit) 
			},silent = TRUE)	
	try({ 
				#allow correlation between beta0 and beta1 but beta2 independent
				tmp.fits[["beta0+beta1+beta2_corf"]] <- tmp.fit <- update( tmp.fit2.e, random=pdBlocked(form=list(beta0+beta1~1,beta2~1),c("pdDiagCorrNeg","pdDiag"))) 
				tmp.aic["beta0+beta2_corf"] <- AIC(tmp.fit) 
			},silent = TRUE)	
}	

.tmp.f <- function(rd.e){
	i_exp <- 31
	i_exp <- "Ul_6_46"
	rde.e <- subset( rd.e, experiment==i_exp )	#entire experiment
	plot(rde.e$resp ~ rde.e$time, pch=19+match(rde.e$replicate,unique(rde.e$replicate)))
	i_rep <- "15"
	rder.e <- subset( rde.e, replicate==i_rep)
	plot(rder.e$resp ~ rder.e$time, pch=19+match(rder.e$replicate,unique(rder.e$replicate)))
	
	tmp.m <- fitKinrespReplicate( rder.e )	
	
	subset(rde.e, resp < 2)
	
	points(rde.e$resp[order(rde.e$time)]~rde.e$time[order(rde.e$time)], col="blue")
	lines(fitted(tmp.fit.e)[order(rde.e$time)]~rde.e$time[order(rde.e$time)], col="red")
}


#---------------- fitting the mixed model to microbial parameters

fitKinrespExpAll <- function(
	### Invoking fitKinrespExp for each experiment within the suite.
	rd.e			##<< dataset with columns suite, experiment, replicate, resp and time, containing only unlimited growth phase.
	, repFits		##<< Initial coefficients mumax, x0, and r0 for each replicate
	, lambda=0.9
	, YCO2=1.5
	, weights=NULL 
){
	# fitKinrespExpAll
	##seealso<< 
	## \code{\link{twKinresp}}
	
	coefExp <- data.frame()
	mExp <- list()
	#i_exp <- "Ek_15_27"
	#i_exp <- "17"
	i_exp <- "3"
	for( i_exp in as.character(unique(rd.e$experiment))){
		print(i_exp)
		rde.e <- subset( rd.e, rd.e$experiment==i_exp )
		tmp.m <- fitKinrespExperiment(rde.e, repFits, lambda=lambda, YCO2=YCO2, weights=weights)$m
		mExp[[ paste(i_exp,sep="_") ]] <- tmp.m
		tmp.coef <- coefKinresp(fixef(tmp.m))
		coefExp <- rbind( data.frame(exp=I(i_exp), as.list(tmp.coef)),coefExp )
	}
	list( m = mExp, coef=coefExp)
}
#mtrace(fitKinrespExpAll)

fitKinrespSuite <- function(
	### Fit microbial form the entire suite.
	rds.e
	, repFits
	, lambda=0.9
	, YCO2=1.5
	, weights=NULL
	, start=NULL
	, tmp.names = c("none","x0l","r0l","x0l+r0l","x0l+r0l+mumaxl")
){
	##seealso<< 
	## \code{\link{twKinresp}}
	
	# XXTodo random effect for each experiment instead of one across all
	# parameters
	#   rds.e: dataset with columns experiment, replicate, resp and time, containing only exponential growht phase
	#   repFits: coefficients mumax, x0, and r0 for each replicate (columns exp)	
	# purpose: estimate the common parameter in the error model 
	# return:
	#	fit: the fit with the lowest AIC may be of class gnls or nlme
	#	aics: informaiton on errors, on various form of inclusion of random effects
	
	if( is.null(start )){
		if( is.null(repFits) ) stop("fitKinrespExpSuite: must provide starting values (start) or estamtes for each replicate (repFits).")
		#tmp.coefs$exp <- factor(tmp.coefs$exp,levels=unique(tmp.coefs$exp))
		#tmp.exp <- as.character(unique(tmp.coefs$exp))
		# take the means of all the replicates as starting values
		#tmp.mumax <- tapply(tmp.coefs[,"mumax"],as.character(tmp.coefs$exp),mean, na.rm=TRUE)
		#tmp.x0l <- log(tapply(tmp.coefs[,"x0"],as.character(tmp.coefs$exp),mean, na.rm=TRUE))
		#tmp.r0l <- logit(tapply(tmp.coefs[,"r0"],as.character(tmp.coefs$exp),mean, na.rm=TRUE))
		#not sure about the order of the starting values. Results of gnls seem to order 10 after 9, so do not use as.character 
		#tmp.mumaxl <- tapply(log(tmp.coefs[,"mumax"]),(tmp.coefs$exp),mean, na.rm=TRUE)
		tmp.exp <- (repFits$experiment)[drop=TRUE] 		
		tmp.mumaxl <- tapply(log(repFits[,"mumax"]),tmp.exp,mean, na.rm=TRUE)
		tmp.x0l <- tapply(log(repFits[,"x0"]),tmp.exp,mean, na.rm=TRUE)
		tmp.r0l <- tapply(logit(repFits[,"r0"]),tmp.exp,mean, na.rm=TRUE)
		start = c( tmp.mumaxl, tmp.r0l, tmp.x0l  ) # r0 comes before x0 in coef
	}
	
	#the same rep might occur in several experiment, but random factor should treat it as differnt, create unqiue rep
	#better use grouped factors 
	#rds.e$exprep = paste( rds.e$experiment, rds.e$replicate, sep="_")
	rds.e$exprep <- with(rds.e, experiment:replicate)[drop=TRUE]
	rds.eg <- groupedData( resp ~ time | experiment/replicate, data=cbind(rds.e, lambda=lambda, YCO2=YCO2))
	# as a parameter tmp.names <- c("none","x0l","r0","x0l+r0","x0l+r0+mumaxl")
	tmp.AIC <- structure( rep(NA, length(tmp.names)),names=tmp.names)
	tmp.fits <- list()
	try({
		tmp.name <- tmp.names[1]
		tmp.fit1 <- gnls(
				resp ~ exp(x0l)*(1-invlogit(r0l))*(1/lambda-1)*exp(mumaxl)/YCO2 + exp(x0l)*invlogit(r0l)*1/lambda*exp(mumaxl)/YCO2 * exp( exp(mumaxl) * time)  
				, params=list(mumaxl+r0l+x0l~experiment)
				, start=start
				, weights=weights
				, data=rds.eg
				, control=gnlsControl(nlsTol=0.1)
				)
		tmp.fits[[tmp.name]] <- tmp.fit1
		tmp.AIC[tmp.name] <- AIC(tmp.fit1)
	})
	#try various models of random effects in different coefficients
	tmp.name <- "x0l"
	for( tmp.name in tmp.names[-1]){
		#tmp.random <- eval(parse(text=paste(tmp.name,"~1|experiment/replicate",sep="")))
		#tmp.random <- eval(parse(text=paste(tmp.name,"~1|exprep",sep="")))
		tmp.random <- eval(parse(text=paste(tmp.name,"~experiment|exprep",sep="")))
		# not accepted as random formula: tmp.random <- eval(parse(text=paste(tmp.name,"~1|experiment:replicate",sep="")))
		#tmp.random <- eval(parse(text=paste("list(experiment=NULL, exprep=",tmp.name,"~1)",sep="")))
		#tmp.random <- eval(parse(text=paste("list(experiment=",1,"~1, exprep=",tmp.name,"~1)",sep="")))
		#tmp.random <- eval(parse(text=paste("list(experiment=NULL, exprep=",tmp.name,"~1)",sep="")))
		print(tmp.random)
		try({
			tmp.fit2 <- nlme(
					resp ~ exp(x0l)*(1-invlogit(r0l))*(1/lambda-1)*exp(mumaxl)/YCO2 + exp(x0l)*invlogit(r0l)*1/lambda*exp(mumaxl)/YCO2 * exp( exp(mumaxl) * time)  
					,fixed=list(mumaxl+r0l+x0l~experiment)	#include fixed effects for each experiment
					,random=tmp.random
					, weights =weights
					, start=start
					, method='REML'	# for unbiased estimates of standard error and confidence intervals 
					,data=rds.eg
					)
			tmp.fits[[tmp.name]] <- tmp.fit2
			tmp.AIC[tmp.name] <- AIC(tmp.fit2)
		})
	}
	tmp.random <- tmp.names[ which.min(tmp.AIC)]
	tmp.fit <- tmp.fits[[ tmp.random  ]]
	list( m=tmp.fit, random=tmp.random, fits=tmp.fits, aics=tmp.AIC )
}

#-------- converting the respiraiton specific outputs to maximum uptake rates and biomass yields
calcKinrespY <- function(
	### Convert yield biomass/respiration to biomass/consumed substrate.	
	YCO2=1.5
){
	##seealso<< 
	## \code{\link{twKinresp}}
	
	(tmp.Y <- 1-1/(1+YCO2))
}

calcKinrespQs <- function(
	### Calculate respiration components.
	mumax
	, YCO2=1.5		##<< Yield per respired CO2
	, lambda=0.9
){
	##seealso<< 
	## \code{\link{twKinresp}}
	
	tmp.Qu = (1/lambda -1)*mumax/YCO2	# maximum uncoupled specific respiration Q' = uncoupled uptake
	tmp.QT = 1/lambda*mumax/YCO2 # Total specific respiration Q_T
	tmp.Y <- calcKinrespY(YCO2)
	tmp.Q = (tmp.QT-tmp.Qu)/(1-tmp.Y)	# maximum substrate uptake for coupled respiration
	list(
		Q=tmp.Q			##<< maximum substrate uptake for coupled respiration
		, Qu=tmp.Qu		##<< maximum uncoupled specific respiration Q' = uncoupled uptake
		, QT=tmp.QT		##<< Total specific respiration Q_T
		, Y=tmp.Y		##<< Yield: biomass/consumed substrate
	)
	### list with components \describe{
	### \item{Q}{maximum substrate uptake for coupled respiration}
	### \item{Qu}{maximum uncoupled specific respiration Q' = uncoupled uptake}
	### \item{QT}{Total specific respiration Q_T}
	### \item{Y}{Yield: biomass/consumed substrate}
	### }	
}



