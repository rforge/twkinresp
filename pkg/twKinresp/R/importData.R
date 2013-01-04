setFactorKinrespData <- function(
	### Codes the columns replicate, experiment, and suite as a factor.
	rd ##<< data.frame holding the experimental observations.
){
	##details<< \describe{\item{Functions for getting the respiration data into correct format}{
	## \itemize{
	## \item{ reading csv-files into a R-variable: \code{\link{read.csv}} \cr Column names must correspond to the variables in \code{\link{assembleKinrespData}} }
	## \item{ coding columns replicate, experiment, and suite as factors: this method }
	## \item{ assembling existing R-variables into correct form: \code{\link{assembleKinrespData}} }
	## \item{ creating a unique identifier for replicates across suits and experiments.: \code{\link{getSERId}} }
	## }
	##}}	

	##seealso<< 
	## \code{\link{twKinresp}}
	
	rd$experiment <- as.factor(rd$experiment)
	rd$replicate <- as.factor(rd$replicate)
	rd$suite <- as.factor(rd$suite)
	rd
	### Modification of input rd.
}
attr(setFactorKinrespData,"ex") <- function(){
    if( FALSE){  # data file not available at R CHECK on R-forge
    	(dataFilename <- file.path( system.file(package="twKinresp"), "data", "respWutzler10.csv" ))
    	# may open the file with a text editor or a spreadsheet program to inspect the format
    	rd0 <- read.csv(dataFilename)	# reading the data from file
    	rd <- setFactorKinrespData(rd0) # Factor coding
    	str(rd)
    	rd[1:10,]
    }
}

assembleKinrespData <- function( 
	### Assemble data into form used by the kinresp package.
	time					##<< Time of respriation measurement
	, resp 					##<< Measured respiration
	, replicate = 1			##<< Identifier of the replicate
	, experiment = ""		##<< Identifier of the experiment: one replicated time series
	, suite= ""				##<< Identifier of the suite: Set of experimental condtions 
){
	##seealso<< 
	## \code{\link{setFactorKinrespData}}
	
	##details<< 
	## All the parameters need to be vectors of the same length.
	## The columens experiment and suite are provided to use mixed models to better constrain
	## those parameters that do not change between experiments.
	data.frame(suite=as.factor(suite), experiment=as.factor(experiment), replicate=as.factor(replicate), time=time, resp=resp)
	### data.frame with colums corresponding to the parameters
}

getSERId <- function(
	### Creating a unique identifier for replicates across suits and experiments.
	rs	##<< the dataset with each row one replicate, with columns suite, experiment, and replicate
){
	paste(rs$suite, rs$experiment, rs$replicate, sep="_")
}






