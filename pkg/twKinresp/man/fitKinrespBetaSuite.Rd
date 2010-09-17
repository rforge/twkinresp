\name{fitKinrespBetaSuite}
\alias{fitKinrespBetaSuite}
\title{fitKinrespBetaSuite}
\description{Fitting the beta-model to all experiments with fixed effects.}
\usage{
fitKinrespBetaSuite(rds.e, repFits, weights = NULL)
}
\arguments{
  \item{rds.e}{dataset with columns experiment, replicate, resp and time, containing only unlimited growth phase (see \code{\link{getUnlimitedGrowthData.kinrespList}}).}
  \item{repFits}{Initial coefficients beta0l,beta1, and beta2l for each replicate}
  \item{weights}{Variance function. see details}
}



\author{Thomas Wutzler <thomas.wutzler@web.de>}





