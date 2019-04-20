\name{kinrespGrowthphaseReplicate}
\alias{kinrespGrowthphaseReplicate}
\title{kinrespGrowthphaseReplicate}
\description{Constrain unlimited growth phase of a single respiration time series.}
\usage{kinrespGrowthphaseReplicate(rder, weights = NULL, 
    orderTime = TRUE, residType = "pearson")}
\arguments{
  \item{rder}{data.frame with columns replicate, time, and resp}
  \item{weights}{Weighting of the observations for non-equal precision}
  \item{orderTime}{if rder is already ordered by time, then you can}
  \item{residType}{type of residuals, see argument type of}
}
\details{It discardsa all data after the inflection point

If the curve deviates from the exponential model, residuals will
be correlated.

The longest time series is selected
for which there is no correlation or a negative correlation
determined by Breusch-Godfrey Test (bgtest) and Durbin-Watson-Test (dwtest)
}
\value{list of class kinresp with components \describe{
\item{dataset}{ the subset with exponential growth phase }
\item{dataGrowth}{ the data of the entire growth phase, i.e. until
maximum respiration rate }
\item{fit}{ the gnls fitting object }
\item{n}{ the number of points suggested by different criteria.
Entry 1 (named "n") gives the best combined estimate.  }
\item{stat}{ the complete statistics r2 and p-values of various
residual tests }
\item{fits}{ results of all the fits. Used e.g. for plotting diagnostics }
}}

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{kinrespGrowthphaseExperiment}}
,\code{\link{plotKinrespDiagnostics.kinresp}}
,\code{\link{twKinresp}}}
\examples{
# we pick and plot the respiration time series of Fig 1 in Wutzler et al. 2010
# data(respWutzler10)
rder <- subset(
  respWutzler10, suite == "Face" & experiment == 3 & replicate == 2 )
plot( resp ~ time, data = rder )

res2 <- kinrespGrowthphaseReplicate(rder, weights = varPower(fixed = 0.5))
res2$n["n"]		#display the number of records
#display the fitting line
lines( fitted(res2$fit) ~ getUnlimitedGrowthData(res2)$time )

# plot diagnostics
plotKinrespDiagnostics(res2)	#use arrow keys to go back in plot history
}
