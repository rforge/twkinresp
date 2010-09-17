\name{coefList.kinrespList}
\alias{coefList.kinrespList}
\title{coefList kinrespList}
\description{Extract the coefficients for each replicate}
\usage{
coefList.kinrespList(tmp.coef, rds.e = NULL, ...)
}
\arguments{
  \item{tmp.coef}{result of \code{\link{kinrespGrowthphaseExperiment}}}
  \item{rds.e}{constrained dataset, which may omit some replicates}
  \item{\dots}{
}
}

\value{named numer matrix (columns experiment, replicate, mumax, x0, and r0) with rows corresponding replicates}

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{coefKinresp.default}}
,\code{\link{twKinresp}}}

