\name{kinrespGrowthphaseExperiment}
\alias{kinrespGrowthphaseExperiment}
\title{kinrespGrowthphaseExperiment}
\description{Confining the unlimited growth phase for several replicates.}
\usage{kinrespGrowthphaseExperiment(rde, ..., stopOnReplicateError = FALSE)}
\arguments{
  \item{rde}{data.frame with columns suite, experiment, replicate, resp and time}
  \item{\dots}{further arguments to \code{\link{kinrespGrowthphaseReplicate}}}
  \item{stopOnReplicateError}{By default replicates where not fit is obtained are ignored.}
}
\details{\describe{\item{Functions for Confining the time series to the unlimited growth phase}{
\itemize{
\item{ Calculating fits and statistics: this method and \code{\link{kinrespGrowthphaseReplicate}} }
\item{ Extracting the records of the unlimited growth phase: \code{\link{getUnlimitedGrowthData.kinresp}} and \code{\link{getUnlimitedGrowthData.kinrespList}} }
\item{ Plotting diagnostic plots: \code{\link{plotKinrespDiagnostics.kinresp}} and \code{\link{plotKinrespDiagnostics.kinrespList}} }
}
}}}
\value{list with components \describe{
\item{n}{data.frame with columns [suite, experiment, replicate], their combination named "ser" and entry n from result of \code{\link{kinrespGrowthphaseReplicate}} for each replicate}
\item{resRep}{list of results of \code{\link{kinrespGrowthphaseReplicate}} for each replicate}
\item{errors}{list of error messages for replicates where fit failed.}
}}

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{twKinresp}}}

