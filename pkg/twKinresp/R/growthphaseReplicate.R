
kinrespGrowthphaseReplicate <- function(
	### Constrain unlimited growth phase of a single respiration time series.
	rder	##<< data.frame with columns replicate, time, and resp
	,weights=NULL		##<< Weighting of the observations for non-equal precision (see details of \code{\link{fitKinrespExperiment}}).
	,orderTime=TRUE		##<< if rder is already ordered by time, then you can set orderTime to FALSE to improve efficiency
	,residType="pearson"	##<< type of residuals, see argument type of \code{\link{residuals.lme}}
){
	##seealso<< 
	## \code{\link{kinrespGrowthphaseExperiment}}
	## ,\code{\link{twKinresp}}
	
	if( length(unique(rder$replicate)) != 1)
		stop("kinrespGrowthphaseReplicate: found other than 1 unique replicate identifier in argument rder.")
	
	# precondition: data is ordered by time
	if( orderTime)
		rder <- rder[order(rder$time),]

	#--------------- discard all points after inflection point ----------
	# got problems with noisy data: include all points up to resp maximum
	#tmp.slope <- diff(c(0,rder$resp))/diff(c(0,rder$time)) 
	#determine the position of the maximum slope (only from 4th to max point 
	#p.infl <- which.max(tmp.slope)
	rder.e <- rder[1:which.max(rder$resp),]
	# formerly also discard also points before the minimum
	#rder.e <- rder.e[which.min(rder.e$resp):nrow(rder.e),]
	# but only discard it analysis of autocorrelation but fit all points in 
	# the beginning else the A parameter will be underestimated or become negative
	#p.min <- which.min(rder.e$resp)
	
	tmp.f <- function(){
		windows(width=4.1,height=3.4,pointsize=10)
		plot(  rder$resp ~ rder$time, type="l" )
		#lot(  tmp.slope ~ rder$time )
		points( rder.e$resp ~ rder.e$time, col="red" )
	}

	#---- fit an exponential curve to subset of less and less points ---
	tmp.n <- floor((nrow(rder.e)/2)) 
	#tmp.n <- 12
	tmp.fits <- list()	#store the model fits
	tmp.bgp <- rep(0,tmp.n)	#p values of the bgtest excluding the points before the minimum
	tmp.bgp1 <- rep(0,tmp.n)	#p values of the bgtest excluding the minimum and the last point
	tmp.bgpfull <- rep(0,tmp.n) #p values of the bgtest 
	tmp.dwp <- rep(0,tmp.n) #p values of the durbin watson dptest, excluding points before minimum 
	tmp.dwp1p <- rep(0,tmp.n) #p values of the durbin watson dptest for positive autocorrelation, excluding pints bef. min and last point
	tmp.dwpfull <- rep(0,tmp.n) #p values of the durbin watson test, all points included
	tmp.dwp1n <- rep(1,tmp.n) #p values of the durbin watson dptest for negative autocorrelation, excluding pints bef. min and last point
	tmp.r2 <- rep(0,tmp.n)
	tmp.r2w <- rep(0,tmp.n)
	tmp.Q <- rep(0,tmp.n)
	i <- 1
	for( i in 1:tmp.n ){
		rdsi = rder.e[1:(nrow(rder.e)+1-i),]
		tmp.fit <- try( fitKinrespBetaReplicate(rdsi$time,rdsi$resp, weights=weights), silent=TRUE )
		if( !inherits(tmp.fit, "try-error")){
			tmp.sd <- try(tmp.fit$sigma*abs(fitted(tmp.fit))^(.twVarFuncCoef(tmp.fit)["power"]))
			tmp.resid = resid(tmp.fit, type=residType)
			#plot( tmp.sd ~ attr(resid(tmp.fit),"std") )
			#plot( tmp.weights ~ rdsi$time)
			#plot( resid(tmp.fit) ~ rdsi$time )
			#plot( resid(tmp.fit, type="pearson") ~ rdsi$time )
			if( inherits(tmp.sd, "try-error")) tmp.sd <- 1
			tmp.fits[[i]] <- tmp.fit
			tmp.bo <- (which.min(rdsi$resp)):(nrow(rdsi))
			tmp.bo1 <- (which.min(rdsi$resp)):(nrow(rdsi)-1)
			tmp.bgpfull[i] <- bgtest(tmp.resid~rdsi$time)$p.value	# in library lmtest
			tmp.bgp[i] <- bgtest(tmp.resid[tmp.bo]~rdsi$time[tmp.bo])$p.value
			tmp.bgp1[i] <- bgtest(tmp.resid[tmp.bo1]~rdsi$time[tmp.bo1])$p.value
			tmp.dwpfull[i] <- dwtest(tmp.resid~rdsi$time, alternative = "greater")$p.value
			tmp.dwp[i] <- dwtest(tmp.resid[tmp.bo]~rdsi$time[tmp.bo], alternative = "greater")$p.value
			tmp.dwp1p[i] <- dwtest(tmp.resid[tmp.bo1]~rdsi$time[tmp.bo1], alternative = "greater")$p.value
			tmp.dwp1n[i] <- dwtest(tmp.resid[tmp.bo1]~rdsi$time[tmp.bo1], alternative = "less")$p.value
			tmp.r2[i] <- 1 - sum(tmp.fit$residuals^2) / sum( (rdsi$resp - mean(rdsi$resp))^2 )  
			tmp.r2w[i] <- 1 - sum(tmp.fit$residuals^2/tmp.sd^2) / sum( (rdsi$resp - mean(rdsi$resp))^2/tmp.sd^2 )  
			tmp.Q[i] <- pchisq( sum(tmp.fit$residuals^2 / tmp.sd^2 ), df=(tmp.fit$dim$N - tmp.fit$dim$p) )   
		}else{
			tmp.fits[[i]] <- NULL
			#all the test statistics are already initialized with autocorrelation p=0 - non food fit
			#dwp1n with 1 to non-significant negative autocorrelation
			#r was initialized with 0: non-good fit
		}
	}
	tmp.fits[[tmp.n+1]] <- "ensure that previous lines are kept"
	
	tmp.stat <-	cbind( n=(nrow(rder.e)+1-(1:tmp.n)), r2=tmp.r2, r2w=tmp.r2w, Q=tmp.Q
			, dwtest1n.p=tmp.dwp1n
			, bgtest1.p=tmp.bgp1, bgtest.p=tmp.bgp, bgtestfull.p=tmp.bgpfull
			#, dwtest1.p=tmp.dwp1p
			, dwtest1=tmp.dwp1p, dwtest.p=tmp.dwp, dwtestfull.p=tmp.dwpfull
			)				
	
	tmp.sl <- 0.05	#significane level
	iseq <- (1:nrow(tmp.stat))
	# determine the i (number of omitted points) for different criteria
	tmp.is <- c( r2 = which.max(tmp.stat[,"r2"]), r2w = which.max(tmp.stat[,"r2w"]), Q = which.max(tmp.stat[,"Q"])
			,dwtest1n.p = {tmp<-which( (tmp.stat[,"dwtest1n.p"][iseq] < tmp.sl)  ); ifelse(length(tmp)>0,min(tmp),Inf) }	#first negative autocorrelation
			,apply(tmp.stat[,-(1:5)],2,function(x){
						i <- {tmp<-which( (x[iseq] >= tmp.sl) & (c(0,x)[iseq] < tmp.sl)); ifelse(length(tmp)>0,min(tmp),Inf)}
					})
			)
	# get the corresponding number of included records
	# suppressWarnings needed for Inf in tmp.is, i.e. when no negative correlation was found.
	tmp.ns <- suppressWarnings( structure( tmp.stat[tmp.is,"n"], names=names(tmp.is) ))	 
	
	i <- min( tmp.is["bgtest1.p"], tmp.is["dwtest1n.p"] )	#cortest
	if( !is.finite(i) ){
		#stop("XXX Warning expPhase: could not determine exponential growth phase: no autocorrelation free fit found.")
		#will give NA in tmp.stat[i,]
	}
	
#	tmp.f <- function(){
#		if( !is.finite(i) )
#			i <- min( tmp.is[c("bgtest.p","bgtestfull.p","dwtest1p.p")] )	#try any other fit by tests
#		#check for the situation where bgest1.p underestimated mu
#		if( !is.na(tmp.ns["bgtest.p"]) && !is.na(tmp.ns["bgtest1.p"]) )
#			if( tmp.ns["bgtest.p"] > tmp.ns["bgtest1.p"]) #usually bgtest1 is less strict so includes more n
#				if( coef(tmp.fits[[tmp.is["bgtest1.p"] ]])[["beta2"]] < coef(tmp.fits[[tmp.is["bgtest.p"] ]])[["beta2"]] )
#					i <- tmp.is["r2"] #see case 29.3 the tests on autocorrelation will not work
#	}	
	
	tmp.ns	<- suppressWarnings(c( cortest=as.numeric(tmp.stat[i,"n"]), tmp.ns))		#append cortest to output of criteria
	
	tmp.ncons <- structure(.kinRespStatN(tmp.ns), names=NULL)	#calc r2wsupp1c
	tmp.ns	<- c( n=tmp.ncons, tmp.ns)					#append r2wsupp1c to output (named n)

	if( is.finite(tmp.ncons) ){
		#rdsi <- rder.e[1:(nrow(rder.e)+1-i),]
		#tmp.fit <- tmp.fits[[i]]
		rdsi <- rder.e[1:tmp.ncons,]
		tmp.fit <- tmp.fits[[ nrow(rder.e)+1-tmp.ncons ]]
	}else{
		#stop("expPhase: could not determine exponential growth phase: no autocorrelation free fit found.")		
		rdsi <- rder.e[FALSE,] # just return empty dataset
		tmp.fit <- NULL
	}
	
	#return
	ret <- list( dataset=rdsi, dataGrowth=rder.e, stat=tmp.stat, fit=tmp.fit, fits=tmp.fits, n=tmp.ns)
	class(ret) <- "kinresp"
	ret
	
	### list of class kinresp with components \describe{
	###  \item{dataset}{ the subset with exponential growth phase }
	###  \item{dataGrowth}{ the data of the entire growth phase, i.e. until maximum respiration rate }
	###  \item{fit}{ the gnls fitting object }
	###  \item{n}{ the number of points suggested by different criteria. Entry 1 (named "n") gives the best combined estimate.  }
	###  \item{stat}{ the complete statistics r2 and p-values of various residual tests }
	###  \item{fits}{ results of all the fits. Used e.g. for plotting diagnostics }
	### }
}
#mtrace(getExpPhase)
#mtrace(getExpPhase,FALSE)

attr(kinrespGrowthphaseReplicate,"ex") <- function(){
	# we pick and plot the respiration time series of Fig 1 in Wutzler et al. 2010
	data(respWutzler10)
	rder <- subset(respWutzler10, suite=="Face" & experiment==3 & replicate==2 )
	plot( resp ~ time, data=rder )
	
	res2 <- kinrespGrowthphaseReplicate(rder, weights=varPower(fixed=0.5)) 
	res2$n["n"]		#display the number of records
	lines( fitted(res2$fit) ~ getUnlimitedGrowthData(res2)$time ) #display the fitting line
}

.kinRespStatN <- function(
	### Combined criterion (r2wsupp1c)  
	tmp.ns		##<< estimated n of basic criteria and cortest
){
	# kinRespStatN
	# here cortest has precedence of r2
	tmp.diffc <- abs(tmp.ns["cortest"] - tmp.ns["r2w"])
	tmp.diffr <- abs(tmp.ns["r2"] - tmp.ns["r2w"])
	# supported by another measure
	if( (!is.na(tmp.diffc) & (tmp.diffc <= 2)) | ( !is.na(tmp.diffr) & (tmp.diffr <= 2))  ){
		#if cortest is near and smaller correct downwards		
		if(  !is.na(tmp.diffc) & (tmp.diffc <= 2) & (tmp.ns["cortest"] < tmp.ns["r2w"])  )
			tmp.ns["r2w"] -1
		else
			tmp.ns["r2w"]
	}else NA
	### adds column n to tmp.ns
}
#kinRespStatN(tmp.ns)


.tmp.f <- function(){
	i_exp <- 11
	rde <- subset(rd, experiment==i_exp)
	i_rep <- 3
	rder <- subset(rde, replicate==i_rep)
	plot(rder$resp ~ rder$time)
	plotFileBasename="output/tmp"
}

