\name{plotKinrespDiagnostics.kinrespList}
\alias{plotKinrespDiagnostics.kinrespList}
\title{plotKinrespDiagnostics kinrespList}
\description{Diagnostic plots for confining unlimited growth phase for each replicate.}
\usage{\method{plotKinrespDiagnostics}{kinrespList}(kinrespRes, plotFileBasename = "", 
    ...)}
\arguments{
  \item{kinrespRes}{object of class kinrespList from \code{\link{kinrespGrowthphaseExperiment}} to plot diagnostics for.}
  \item{plotFileBasename}{basename of figure files for diagnostic plots}
  \item{\dots}{further argument to \code{\link{plotKinrespDiagnostics.kinresp}}}
}

\value{invisible: list of all results of \code{\link{plotKinrespDiagnostics.kinresp}} }

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{kinrespGrowthphaseExperiment}}
,\code{\link{plotKinrespDiagnostics.kinresp}}
,\code{\link{twKinresp}}}
\examples{FALSE}
