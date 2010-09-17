\name{plotKinrespDiagnostics.kinresp}
\alias{plotKinrespDiagnostics.kinresp}
\title{plotKinrespDiagnostics kinresp}
\description{Diagnostic plots for confining unlimited growth phase.}
\usage{
plotKinrespDiagnostics.kinresp(kinrespRes, plotFileBasename = "", expText = NULL, xunit = "hour", yunit = "gC CO2 / gC soil / hour", kind = c(resp = "resp", respLog = "respLog", resid = "resid"), residType = "pearson", rangeN = nrow(kinrespRes$dataGrowth) + c(1 - nrow(kinrespRes$stat), 0), ...)
}
\arguments{
  \item{kinrespRes}{object of class kinresp from \code{\link{kinrespGrowthphaseReplicate}} to plot diagnostics for.}
  \item{plotFileBasename}{basename of figure files for diagnostic plots}
  \item{expText}{
}
  \item{xunit}{unit of time to be displayed in the x-axis title}
  \item{yunit}{unit of respiration to be displayed in the y-axis title}
  \item{kind}{kind of diagnostic plot series: character vector of resp, respLog, and resid\describe{
\item{resp}{Observations and fitted model at original scale.}
\item{respLog}{Observations and fitted model at logarithmic scale.}
\item{resid}{Model residuals}
}}
  \item{residType}{type of residuals, see argument type of \code{\link{residuals.lme}}}
  \item{rangeN}{the range (numeric vector of length two) of the number of records, for which diagnostic plots should be produced.}
  \item{\dots}{
}
}
\details{\describe{\item{Writing figure files}{
If plotFileBasename has been specified, diagnostic plot series are save in   
emf and pdf format
}}}
\value{invisible(NULL)}

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{kinrespGrowthphaseExperiment}}
,\code{\link{plotKinrespDiagnostics.kinrespList}}
,\code{\link{twKinresp}}}

