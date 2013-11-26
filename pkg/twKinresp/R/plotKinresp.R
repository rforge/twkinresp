setMethodS3("plotKinrespDiagnostics","default", function( 
	### Diagnostic plots for confining unlimited growth phase.
		kinrespRes
		,...
	){
		##seealso<< 
		## ,\code{\link{plotKinrespDiagnostics.kinresp}}
		## ,\code{\link{twKinresp}}
		
		stop("plotKinrespDiagnostics.default: unknown class of argument kinrespRes. Must be a result of kinrespGrowthphaseReplicate or kinrespGrowthphaseExperiment")
	})

setMethodS3("plotKinrespDiagnostics","kinresp", function( 
	### Diagnostic plots for confining unlimited growth phase.
	kinrespRes				##<< object of class kinresp from \code{\link{kinrespGrowthphaseReplicate}} to plot diagnostics for.
	,plotFileBasename=""	##<< basename of figure files for diagnostic plots
	,expText=NULL
	,xunit="hour"			##<< unit of time to be displayed in the x-axis title
	,yunit="gC CO2 / gC soil / hour"	##<< unit of respiration to be displayed in the y-axis title
	,kind= 		##<< kind of diagnostic plot series: character vector of resp, respLog, and resid 
		##describe<<  
		c(resp="resp"	   ##<< Observations and fitted model at original scale.
		,respLog="respLog" ##<< Observations and fitted model at logarithmic scale.
		,resid="resid"	   ##<< Model residuals
		)
		##end<< 
	,residType="pearson"   ##<< type of residuals, see argument type of \code{\link{residuals.lme}}
	,rangeN=nrow(kinrespRes$dataGrowth)+c(1-nrow(kinrespRes$stat),0)
		### the range (numeric vector of length two) of the number of records, for which diagnostic plots should be produced.
	,...	##<< further  arguments to array plotting 
		## \itemize{
		## \item rangeN: nCol integer scalar of number of columns
		## \item topOffsetLine=0.08: inset of the diagnostics 
		## }
){
	# plotKinrespDiagnostics.kinresp
		##seealso<< 
		## \code{\link{kinrespGrowthphaseReplicate}}
		## ,\code{\link{kinrespGrowthphaseExperiment}}
		## ,\code{\link{plotKinrespDiagnostics.kinrespList}}
		## ,\code{\link{twKinresp}}
		
	##details<< \describe{\item{Writing figure files}{
	## If plotFileBasename has been specified, diagnostic plot series are save in   
	## emf and pdf format
	## }}
	
	rder.e <- kinrespRes$dataGrowth
	#constrain rangeN to possible values
	tmp.rangeNposs <- nrow(kinrespRes$dataGrowth)+c(1-nrow(kinrespRes$stat),0)
	tmp.rangeNArg <- sort(rangeN[1:2])	
	tmp.rangeN <- c( max(tmp.rangeNposs[1],tmp.rangeNArg[1]), min(tmp.rangeNposs[2],tmp.rangeNArg[2]) )
	tmp.fits <- kinrespRes$fits
	tmp.stat <- kinrespRes$stat
	tmp.is <- nrow(rder.e)+1 - kinrespRes$n
	#mtrace(plotRepResidualArray)
	if( "resp" %in% kind) .plotRepRespArray (rder.e,tmp.fits,tmp.stat,tmp.is,plotFileBasename=plotFileBasename,expText=expText, xunit=xunit, yunit=yunit, tmp.rangeN, ...)
	if( "respLog" %in% kind) .plotRepRespLogArray (rder.e,tmp.fits,tmp.stat,tmp.is,plotFileBasename=plotFileBasename,expText=expText, xunit=xunit, yunit=yunit, tmp.rangeN, ...)
	if( "resid" %in% kind) .plotRepResidualArray (rder.e,tmp.fits,tmp.stat,tmp.is,plotFileBasename=plotFileBasename,expText=expText, xunit=xunit, yunit=yunit, residType=residType, tmp.rangeN, ...)
	return(invisible(NULL))
	### invisible(NULL)
})

setMethodS3("plotKinrespDiagnostics","kinrespList", function(
	### Diagnostic plots for confining unlimited growth phase for each replicate.
	kinrespRes				##<< object of class kinrespList from \code{\link{kinrespGrowthphaseExperiment}} to plot diagnostics for.
	,plotFileBasename=""	##<< basename of figure files for diagnostic plots
	,...					##<< further argument to \code{\link{plotKinrespDiagnostics.kinresp}}
){
	# plotKinrespDiagnostics.kinrespList
	##seealso<< 
	## \code{\link{kinrespGrowthphaseExperiment}}
	## ,\code{\link{plotKinrespDiagnostics.kinresp}}
	## ,\code{\link{twKinresp}}
	res <- lapply(kinrespRes$resRep,function(kinrespRep){
		plotFileBasenameRep <- ifelse( length(plotFileBasename) > 0
			,paste(plotFileBasename,kinrespRep$dataset$experiment[1],kinrespRep$dataset$replicate[1],sep="_")
			,"")
		plotKinrespDiagnostics.kinresp(kinrespRep,plotFileBasename=plotFileBasenameRep,...)
	}) 
	return(invisible(res))
	### invisible: list of all results of \code{\link{plotKinrespDiagnostics.kinresp}} 
})

.plotRepResidualArray <- function(rder.e,tmp.fits,tmp.stat,tmp.is,plotFileBasename="",expText=NULL, residType="pearson", xunit="hour", yunit="gC CO2 / gC soil / hour"
	,rangeN
	,nCol=NA
	,topOffsetLine=0.08
){
	# rder.e: the data of the one replicate that was used for i=1 (including all points)
	# tmp.n: the number of fits: 1 all points 2 one point omitted ..
	# tmp.fits: the actual fitting objects with residuals 
	# tmp.stat: test statistics
	# tmp.is: number of winning fit by criterion
	# plotFileBasename: the basename of the output file
	# expTest: additional test printed along the x axis, usually denoting the experiment and replicate
	# if the model was gls or nlme model, then the weighted residuals are plotted (see residuals.gls)
	#windows(pointsize=10,width=6.4,height=6.4)
	tmp.n <- diff(rangeN)+1
	tmp <- if(is.finite(nCol)) nCol else ceiling(sqrt(tmp.n))
	par( 
		#mfrow=c(4,3)
		mfrow=c(ceiling(tmp.n/tmp), tmp)
		, mar=c(0, 0, 0, 0) + 0.1 
		, oma=c(1.5, 1.5, 0, 0) + 0.1 
	)
	nmax <- nrow(rder.e)
	for( i in nmax+1-(rangeN[2]:rangeN[1]) ){
		rdsi = rder.e[1:(nmax+1-i),]
		tmp.fit <- tmp.fits[[i]]
		if( is.null(tmp.fit) ){
			plot(1~1,type="n",axes=FALSE, xlab="", ylab="")
			box()
		}else{
			tmp.r <- resid( tmp.fit, type=residType ) #if it is a gls fits, the weighted (standardized) 
			plot( tmp.r ~ rdsi$time, axes=FALSE, col="grey70" )
			box()
			mtext(paste(
					"n=",nrow(rder.e)+1-i
					,sep=""),side=3,line=-1.1-topOffsetLine, adj=0.02)
			tmp.c <- "bgtestfull.p"
			mtext(paste("p.1=",round(tmp.stat[i,tmp.c],3)
					,sep=""),side=1,line=-3.3-topOffsetLine, adj=0.02
				, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
			tmp.c <- "bgtest.p"
			mtext(paste("p.2=",round(tmp.stat[i,tmp.c],3)
					,sep=""),side=1,line=-2.2-topOffsetLine, adj=0.02
				, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
			tmp.c <- "bgtest1.p"
			mtext(paste("p.3=",round(tmp.stat[i,tmp.c],3)
					,sep=""),side=1,line=-1.1-topOffsetLine, adj=0.02
				, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
			tmp.c <- "dwtest1n.p"
			mtext(paste("p.n=",round(tmp.stat[i,"dwtest1n.p"],3)
					,sep=""),side=1,line=-1.1-topOffsetLine, adj=0.98
				, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
		}
	}
	if( is.null(expText) ) 
		expText <- paste( "(Experiment ",rder.e$experiment[1], "_",rder.e$replicate[1], ")", sep="" )
	tmp.time <- paste(format(range(rder.e$time),digits=2,trim=TRUE),collapse=" to ")
	mtext( paste("Time",tmp.time, xunit,expText),1, line=0.3,outer=TRUE)
	#tmp.resp <- paste(format(range(tmp.r),digits=2,trim=TRUE),collapse=" to ")
	mtext( paste("Respiration residual", yunit),2,line=0.3,outer=TRUE,las=0)
	if( plotFileBasename != ""){
		tmp.t<-paste(plotFileBasename,"_residArray",sep="");savePlot(tmp.t,type="emf");savePlot(tmp.t,type="eps");savePlot(tmp.t,type="pdf")
	}
}

.plotRepRespArray <- function(rder.e,tmp.fits,tmp.stat,tmp.is,plotFileBasename="",expText=NULL, xunit="hour", yunit="gC CO2 / gC soil / hour"
	,rangeN
	,nCol=NA
	,topOffsetLine=0.08
){
	tmp.n <- diff(rangeN)+1
	# rder.e: the data of the one replicate that was used for i=1 (including all points)
	# tmp.n: the number of fits: 1 all points 2 one point omitted ..
	# tmp.fits: the actual fitting objects with residuals 
	# tmp.stat: test statistics
	# tmp.is: number of winning fit by criterion
	# plotFileBasename: the basename of the output file
	# expTest: additional test printed along the x axis, usually denoting the experiment and replicate
	
	#windows(pointsize=10,width=6.4,height=6.4)
	tmp <- if(is.finite(nCol)) nCol else ceiling(sqrt(tmp.n))
	par( 
		#mfrow=c(4,3)
		mfrow=c(ceiling(tmp.n/tmp), tmp)
		, mar=c(0, 0, 0, 0) + 0.1 
		, oma=c(1.5, 1.5, 0, 0) + 0.1 
	)
	nmax <- nrow(rder.e)
	for( i in nmax+1-(rangeN[2]:rangeN[1]) ){
		rdsi = rder.e[1:(nmax+1-i),]
		tmp.fit <- tmp.fits[[i]]
		
		plot( rder.e$resp ~ rder.e$time, axes=FALSE, col="grey70" )
		box()
		points( rdsi$resp ~ rdsi$time )
		if( !is.null(tmp.fit) ){
			#extract the confidence interval of beta2 (mumax)
			tmp2 <- .confintKinrespBeta(confint(tmp.fit))["beta2",]
			#tmp <- summary(tmp.fit)$tTable["beta2",2]
			#tmp2 <- coef(tmp.fit)[["beta2"]]+1.96*c(-tmp,+tmp)
			tmp3 <- signif(tmp2,3)
			tmp4 <- paste("mu=(", tmp3[1],",",tmp3[2],")",sep="")
			
			lines( predict(tmp.fit,newdata=data.frame(x=rder.e$time,y=rder.e$resp)) ~ rder.e$time)
			mtext(paste("n=",nrow(rder.e)+1-i,sep=""),side=3,line=-1.1-topOffsetLine, adj=0.02)
			# lines -1.1, -2.3, -3.5, -4.7,  -5.9
			mtext(tmp4,side=3,line=-2.3-topOffsetLine, adj=0.02)
			tmp.c <- "r2"
			mtext(paste("r2=",signif(tmp.stat[i,tmp.c],5),sep="")
				,side=3,line=-3.5-topOffsetLine, adj=0.02
				, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
			tmp.c <- "r2w"
			mtext(paste("r2w=",signif(tmp.stat[i,tmp.c],5),sep="")
				,side=3,line=-4.7-topOffsetLine, adj=0.02
				, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
			tmp.c <- "Q"
			mtext(paste("Q=",signif(tmp.stat[i,tmp.c],5),sep="")
				,side=3,line=-5.9-topOffsetLine, adj=0.02
				, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
			#mtext(paste("delta=",signif(.twVarFuncCoef(tmp.fit),3),sep=""),side=3,line=-5.9, adj=0.02)
			#mtext(paste("mu=",signif(coef(tmp.fit)[["beta2"]],3),sep=""),side=3,line=-5.9, adj=0.02)
		}
	}
	if( is.null(expText) ) 
		expText <- paste( "(Experiment ",rder.e$experiment[1], "_",rder.e$replicate[1], ")", sep="" )
	tmp.time <- paste(format(range(rder.e$time),digits=2,trim=TRUE),collapse=" to ")
	mtext( paste("Time",tmp.time, xunit,expText),1, line=0.3,outer=TRUE)
	tmp.resp <- paste(format(range(rder.e$resp),digits=2,trim=TRUE),collapse=" to ")
	mtext( paste("Respiration Rate",tmp.resp, yunit),2,line=0.3,outer=TRUE,las=0)
	if( plotFileBasename != ""){
		tmp.t<-paste(plotFileBasename,"_respArray",sep="");savePlot(tmp.t,type="emf");savePlot(tmp.t,type="eps");savePlot(tmp.t,type="pdf")
	}
}
#mtrace(plotRepRespArray)

.plotRepRespLogArray <- function(rder.e,tmp.fits,tmp.stat,tmp.is,plotFileBasename="",expText=NULL, xunit="hour", yunit="gC CO2 / gC soil / hour"
	,rangeN
	,nCol=NA
	,topOffsetLine=0.08
){
	tmp.n <- diff(rangeN)+1
	# rder.e: the data of the one replicate that was used for i=1 (including all points)
	# tmp.n: the number of fits: 1 all points 2 one point omitted ..
	# tmp.fits: the actual fitting objects with residuals 
	# tmp.stat: test statistics
	# tmp.is: number of winning fit by criterion
	# plotFileBasename: the basename of the output file
	# expTest: additional test printed along the x axis, usually denoting the experiment and replicate
	
	#windows(pointsize=10,width=6.4,height=6.4)
	tmp <- if(is.finite(nCol)) nCol else ceiling(sqrt(tmp.n))
	par( 
		#mfrow=c(4,3)
		mfrow=c(ceiling(tmp.n/tmp), tmp)
		, mar=c(0, 0, 0, 0) + 0.1 
		, oma=c(1.5, 1.5, 0, 0) + 0.1 
	)
	# determine beta0 for each constrained dataset
	#sapply( 1:(length(tmp.fits)-1), function(i){ print(i); if(!is.null(tmp.fits[[i]])) coefKinrespBeta(coef(tmp.fits[[i]]))[["beta0"]] else NA})
	#tmp.fit <- tmp.fits[[6]]]
	beta0 <- sapply( 1:(length(tmp.fits)-1), function(i){ if(!is.null(tmp.fits[[i]])) coefKinrespBeta(coef(tmp.fits[[i]]))[["beta0"]] else NA})  
	beta0All <- min(beta0, na.rm=TRUE)	#lenght-1 because last is only placeholder
	beta0[ is.na(beta0)] <- beta0All
	
	#coefKinresp(coef(tmp.fits[[16]]))[["x0"]]
	nmax <- nrow(rder.e)
	for( i in nmax+1-(rangeN[2]:rangeN[1]) ){
		rdsi = rder.e[1:(nmax+1-i),]
		tmp.fit <- tmp.fits[[i]]
		tmp.logDiff <- suppressWarnings(log(rder.e$resp-beta0[i]))
		plot( tmp.logDiff ~ rder.e$time, axes=FALSE, col="grey70" ) 	
		box()
		tmp.logDiffi <- suppressWarnings(log(rdsi$resp-beta0[i]))
		points( tmp.logDiffi ~ rdsi$time )
		if( !is.null(tmp.fit) ){
			#extract the confidence interval of beta2 (mumax)
			tmp2 <- .confintKinrespBeta(confint(tmp.fit))["beta2",]
			#tmp <- summary(tmp.fit)$tTable["beta2",2]
			#tmp2 <- coef(tmp.fit)[["beta2"]]+1.96*c(-tmp,+tmp)
			tmp3 <- signif(tmp2,3)
			tmp4 <- paste("mu=(", tmp3[1],",",tmp3[2],")",sep="")
			
			lines( log(predict(tmp.fit,newdata=data.frame(x=rder.e$time,y=rder.e$resp))-beta0[i]) ~ rder.e$time)
			mtext(paste("n=",nrow(rder.e)+1-i,sep=""),side=3,line=-1.1-topOffsetLine, adj=0.02)
			# lines -1.1, -2.3, -3.5, -4.7,  -5.9
			mtext(tmp4,side=3,line=-2.3, adj=0.02)
			#tmp.c <- "r2"
			#mtext(paste("r2=",signif(tmp.stat[i,tmp.c],5),sep="")
		#		,side=3,line=-3.5-topOffsetLine, adj=0.02
	#			, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
			tmp.c <- "r2w"
			mtext(paste("r2w=",signif(tmp.stat[i,tmp.c],4),sep="")
				,side=3,line=-3.5-topOffsetLine, adj=0.02
				, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
			tmp.c <- "Q"
			mtext(paste("Q=",signif(tmp.stat[i,tmp.c],5),sep="")
				,side=3,line=-4.7-topOffsetLine, adj=0.02
				, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
			#tmp.c <- "r2w"
			#mtext(paste("r2w=",signif(tmp.stat[i,tmp.c],5),sep="")
			#			,side=3,line=-3.5, adj=0.02
			#				, font={tmp <- 1; if(!is.na(tmp.is[tmp.c]) & tmp.is[tmp.c] == i){ tmp <- 2}; tmp})
			#mtext(paste("delta=",signif(.twVarFuncCoef(tmp.fit),3),sep=""),side=3,line=-5.9, adj=0.02)
			#mtext(paste("mu=",signif(coef(tmp.fit)[["beta2"]],3),sep=""),side=3,line=-5.9, adj=0.02)
		}
	}
	if( is.null(expText) ) 
		expText <- paste( "(Experiment ",rder.e$experiment[1], "_",rder.e$replicate[1], ")", sep="" )
	tmp.time <- paste(format(range(rder.e$time),digits=2,trim=TRUE),collapse=" to ")
	mtext( paste("Time",tmp.time, xunit,expText),1, line=0.3,outer=TRUE)
	tmp.resp <- paste(format(range(rder.e$resp),digits=2,trim=TRUE),collapse=" to ")
	mtext( paste("log(Respiration Rate",tmp.resp, yunit,"- beta0)"),2,line=0.3,outer=TRUE,las=0)
	if( plotFileBasename != ""){
		tmp.t<-paste(plotFileBasename,"_respLogArray",sep="");savePlot(tmp.t,type="emf");savePlot(tmp.t,type="eps");savePlot(tmp.t,type="pdf")
	}
}
#mtrace(plotRepRespLogArray)

