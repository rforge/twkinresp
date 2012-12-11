\name{setFactorKinrespData}
\alias{setFactorKinrespData}
\title{setFactorKinrespData}
\description{Codes the columns replicate, experiment, and suite as a factor.}
\usage{setFactorKinrespData(rd)}
\arguments{
  \item{rd}{data.frame holding the experimental observations.}
}
\details{\describe{\item{Functions for getting the respiration data into correct format}{
\itemize{
\item{ reading csv-files into a R-variable: \code{\link{read.csv}} \cr Column names must correspond to the variables in \code{\link{assembleKinrespData}} }
\item{ coding columns replicate, experiment, and suite as factors: this method }
\item{ assembling existing R-variables into correct form: \code{\link{assembleKinrespData}} }
\item{ creating a unique identifier for replicates across suits and experiments.: \code{\link{getSERId}} }
}
}}	}
\value{Modification of input rd.}

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{twKinresp}}}
\examples{
(dataFilename <- file.path( system.file(package="twKinresp"), "data", "respWutzler10.csv" ))
# may open the file with a text editor or a spreadsheet program to inspect the format
rd0 <- read.csv(dataFilename)	# reading the data from file
rd <- setFactorKinrespData(rd0) # Factor coding
str(rd)
rd[1:10,]
}
