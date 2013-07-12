#--------------------- growth phase of experiment ---------------------
kinrespGrowthphaseExperiment <- function(
	### Confining the unlimited growth phase for several replicates.
	rde		##<< data.frame with columns suite, experiment, replicate, resp and time
	, ...	##<< further arguments to \code{\link{kinrespGrowthphaseReplicate}}
	, stopOnReplicateError=FALSE	##<< By default replicates where not fit is obtained are ignored.
){
	##seealso<< 
	## \code{\link{twKinresp}}
	
	##details<< \describe{\item{Functions for Confining the time series to the unlimited growth phase}{
	## \itemize{
	## \item{ Calculating fits and statistics: this method and \code{\link{kinrespGrowthphaseReplicate}} }
	## \item{ Extracting the records of the unlimited growth phase: \code{\link{getUnlimitedGrowthData.kinresp}} and \code{\link{getUnlimitedGrowthData.kinrespList}} }
	## \item{ Plotting diagnostic plots: \code{\link{plotKinrespDiagnostics.kinresp}} and \code{\link{plotKinrespDiagnostics.kinrespList}} }
	## }
	##}}	
	
	if( length(unique(rde$suite)) != 1)
		stop("kinrespGrowthphaseExperiment: found other than 1 unique suite identifier in argument rde")
	if( length(unique(rde$experiment)) != 1)
		stop("kinrespGrowthphaseExperiment: found other than 1 unique experiment identifier in argument rde")
	
	repIds <- as.character(unique(rde$replicate))
	expName <- rde$experiment[1]
	suiteName <- rde$suite[1]
	exprepIds <- paste(as.character(suiteName),as.character(expName),as.character(repIds),sep="_")
	tmp.errs <- list()
	tmp.n <- list()	#will be converted to data.frame afterwards
	tmp.resRep <- list()
	
	i <- 2
	#get the experimental growth phase of each replicate
	for( i in seq_along(repIds)){
		i_exprep = exprepIds[i]
		i_rep = repIds[i]
		print(paste("ser(suite_experiment_replicate)=",i_exprep,sep=""))
		rder <- subset(rde, replicate==i_rep)
		tmp.res <- if( stopOnReplicateError ){
			tmp <- kinrespGrowthphaseReplicate(rder, ... )
		}else{
			tmp.res <- try({ # this time with try
					tmp <- kinrespGrowthphaseReplicate(rder, ...	)
			})
		}
		# record an eventual error, but go on with the other replicates
		if( inherits(tmp.res, "try-error")){
			tmp.errs[[i_exprep]] <- tmp.res[1]
		}else{
			tmp.resRep[[i_exprep]] <- tmp.res 
			tmp.n[[i_exprep]] <- tmp.res$n
		}
	}
	
	# if no replicate was fitted (no data) signal an error
	#   dataset that only contains the exponential growth phase of each replicate
	#	errors: informaiton on errors, on some replicates
	dfn <- do.call( rbind, tmp.n )
	dfn2 <- cbind( data.frame( suite=suiteName, experiment=expName, replicate=repIds, ser=names(tmp.resRep)), dfn )
	res <- list(
		n=dfn2				
		, resRep=tmp.resRep	 
		, errors=tmp.errs	 
	)
	class(res) <- "kinrespList"
	res
	### list with components \describe{
	### \item{n}{data.frame with columns [suite, experiment, replicate], their combination named "ser" and entry n from result of \code{\link{kinrespGrowthphaseReplicate}} for each replicate}
	### \item{resRep}{list of results of \code{\link{kinrespGrowthphaseReplicate}} for each replicate}
	### \item{errors}{list of error messages for replicates where fit failed.}
	### }
}
#mtrace(kinrespGrowthphaseExperiment)

combineKinrespLists <- function( 
	### Combine several kinrespList to one bigger kinrespList.
	x		##<< List of kinrespList entries
){
	if( 0==length(x)) return( NULL )
	if( !is.list(x) ) {warning("combineKinrespLists argument must be a list of kinrespList entries, Returning NULL"); return(NULL) }
	if( !is(x[[1]],"kinrespList") ) {warning("combineKinrespLists argument must be a list of kinrespList entries, Returning NULL"); return(NULL) }
	# need to remove the names on lapply, else the indices will prefix the names of the entries
	res <- list(
		n = do.call( rbind, lapply(x,"[[","n") )
		,resRep = tmp <- do.call( c,  {tmp<-lapply(x,"[[","resRep");names(tmp)<-NULL;tmp} )
		,errors = do.call( c,  {tmp<-lapply(x,"[[","errors");names(tmp)<-NULL;tmp} )
		)
	rownames(res$n) <- NULL
	class(res) <- "kinrespList"
	res
}

R.methodsS3::setMethodS3("getUnlimitedGrowthData","default", function( 
	### Extract the dataset with records of unlimited growth phase.
	kinrespRes
	,...
){
	##seealso<< 
	## \code{\link{getUnlimitedGrowthData.kinresp}}
	## \code{\link{getUnlimitedGrowthData.kinrespList}}
	stop("getUnlimitedGrowthData: unknown class of argument kinrespRes. Must be a result of kinrespReplicate or kinrespGrowthphaseExperiment")
})
R.methodsS3::setMethodS3("getUnlimitedGrowthData","kinresp", function( 
	### Extract the dataset with records of unlimited growth phase.
	kinrespRes				##<< object of class kinresp from \code{\link{kinrespGrowthphaseReplicate}}.
	,n=kinrespRes$n["n"]	##<< the number of records in the growth phase
	,...
){
	# getUnlimitedGrowthData.kinresp
	##seealso<< 
	## \code{\link{kinrespGrowthphaseExperiment}}
	## ,\code{\link{getUnlimitedGrowthData.kinrespList}}
	## ,\code{\link{twKinresp}}
		
	kinrespRes$dataGrowth[1:n,]
})

R.methodsS3::setMethodS3("getUnlimitedGrowthData","kinrespList", function(
	### Extract the dataset with records of unlimited growth phase.
	kinrespRes	    ##<< object of class kinrespList from \code{\link{kinrespGrowthphaseExperiment}}.
	,n=integer(0)	##<< integer vector with names (character suite_experiment_replicate) speciying the number of records in growth phase, 
        ##<< if not given (length zero), the defaults from \code{kinrespRes$n["n"]} are used 
	,...
){
	# getUnlimitedGrowthData.kinrespList
	##seealso<< 
	## \code{\link{kinrespGrowthphaseExperiment}}
	## ,\code{\link{getUnlimitedGrowthData.kinresp}}
	## ,\code{\link{twKinresp}}
	
	# kinrespRepName <- names(kinrespRes$resRep)[1]
	# n[kinrespRepName] = 12
    if( length(n) ){ # number are explicitly given
        if( !is.numeric(n) ) stop("getUnlimitedGrowthData.kinrespList: n must be an integer vector.")
        if( length(names) ){
            if( !all(names(n) == names(kinrespRes$resRep)) ) stop("getUnlimitedGrowthData.kinrespList: n has wrong names. Either give no names or give names corresponding to a name of the replicate series.")    
        }else
            if( length(n) != length(kinrespRes$resRep) ) stop("getUnlimitedGrowthData.kinrespList: n has wrong length. It must have an entry for each replicate.")
    }else{
        # defaults from replicates
        n <- sapply( kinrespRes$resRep, function(kinrespRep){ kinrespRep$n["n"]} )
        names(n) <- names(kinrespRes$resRep)
    }
	#kinrespRepName <- names(kinrespRes$resRep)[1]
	res <- do.call( rbind, lapply(names(kinrespRes$resRep),function(kinrespRepName){
		kinrespRep <- kinrespRes$resRep[[kinrespRepName]]
        nRep <- n[kinrespRepName]
		if( is.finite(nRep) )		
			getUnlimitedGrowthData( kinrespRep, n=nRep )
		else
			kinrespRep$dataGrowth[FALSE,]
	}) )
	rownames(res) <- NULL
	res
})


#------------ Application for each experiment in the dataset
kinrespGrowthphase <- function(
	### Confining the unlimited growth phase all replicates in the dataset.
	rd		##<< data.frame with columns suite, experiment, replicate, resp and time
	, ...	##<< further arguments to \code{\link{kinrespGrowthphaseExperiment}}
){
	##seealso<< 
	## \code{\link{twKinresp}}
	# experiment identifier may be the same for different suits
	suiteExp <- (rd$suite:rd$experiment)[drop=TRUE]	#dropping unused levels see ?interaction
	#rde <- subset(rd,suiteExp=="Face:3") 
	resl <- by( rd, suiteExp, function(rde){
		#resExp <- resExp0 <- kinrespGrowthphaseExperiment(rde,weights=varPower(fixed=0.3))
		#kinrespGrowthphaseExperiment(rde,weights=varPower(fixed=0.3),stopOnReplicateError=TRUE)		
		kinrespGrowthphaseExperiment(rde,...)		
	})
	names(resl) <- rep("", length(resl))
	res <- list(
		n = do.call( rbind, lapply( resl, function(resExp){resExp[["n"]]}) )
		,resRep = do.call( c, lapply( resl, function(resExp){resExp[["resRep"]]}) )
		,errors = do.call( c, lapply( resl, function(resExp){resExp[["errors"]]}) )
	)
	class(res) <- class(resl[[1]])
	res
	### same format as \code{\link{kinrespGrowthphaseExperiment}}
}
#mtrace(kinrespGrowthphase)



