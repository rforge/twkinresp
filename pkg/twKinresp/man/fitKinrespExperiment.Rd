\name{fitKinrespExperiment}
\alias{fitKinrespExperiment}
\title{fitKinrespExperiment}
\description{Fit microbial form to several replicated respiration time series of one experiment.}
\usage{fitKinrespExperiment(rde.e, repFits, lambda = 0.9, YCO2 = 1.5, 
    weights = NULL, start = NULL, tmp.names = c("none", "x0l", 
        "r0l", "x0l+r0l", "x0l+r0l+mumaxl"), showFitErrorMsg = FALSE)}
\arguments{
  \item{rde.e}{dataset with columns experiment, replicate, resp and time, containing only unlimited growth phase (see \code{\link{getUnlimitedGrowthData.kinrespList}}).}
  \item{repFits}{Initial coefficients mumax, x0, and r0 for each replicate (see \code{\link{coefKinresp.kinrespList}}).}
  \item{lambda}{Ratio of growth associated (coupled) specific respiration to total specific respiration. Usually 0.9.}
  \item{YCO2}{Ratio of assimilated carbon per respired carbon (Y/(1-Y)).  Usually 1.5.}
  \item{weights}{Variance function. see details}
  \item{start}{Specifying starting values in an alternative way to argument repFits.}
  \item{tmp.names}{scenarios of random effects}
  \item{showFitErrorMsg}{if FALSE (standard) then error msg are suppressed when no fit is obtained. This is a common case for the test of variants that include too many random effects.}
}
\details{\describe{\item{Microbial parameters across replicates}{
Microbial parameters can be inferred from fitting against the time series of respiration 
of each replicate separately (\code{\link{fitKinrespReplicate}}.
However, what are then parameters of the population? The average across replicate parameters will be wrong.

Alternatively, one can first average the measurement across replicates at for each measurement time.
However, the uncertainty of the parameters will be wrong.

A viable soluion is to fit a mixed model to all the replicate data. The micoribal parameters, then are
described by a mean value of the population and a variance across replicates.
Several options of which parameters vary across replicates can be tested.
}}

\describe{\item{Variance function}{
If no weights are given, the measurement errors of the single
observations are assumed to be identical. If measurement errors
increase with the magnitude of the observations, this can be
modeled by power variance function:
\code{weights=\link{varPower}(fixed=delta)},
with delta being a value between 0 (constant expected absolute error) 
and 1 (constant expected relative error).
}}}
\value{list with components \describe{
\item{model}{best fitting result}
\item{random}{best random effects scenario}
\item{fits}{fits of all random effects scenarios}
\item{aics}{Akaike information criterion for random effects scenarios}
}}

\author{Thomas Wutzler <thomas.wutzler@web.de>}




\examples{
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
windows(record=TRUE)
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
