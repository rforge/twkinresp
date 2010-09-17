\name{getUnlimitedGrowthData.kinrespList}
\alias{getUnlimitedGrowthData.kinrespList}
\title{getUnlimitedGrowthData kinrespList}
\description{Extract the dataset with records of unlimited growth phase.}
\usage{
getUnlimitedGrowthData.kinrespList(kinrespRes, n = list(), ...)
}
\arguments{
  \item{kinrespRes}{object of class kinrespList from \code{\link{kinrespGrowthphaseExperiment}}.}
  \item{n}{named integer vector with names (character suite_experiment_replicate) speciying the number of records in growth phase, overriding the one obtained from \code{kinrespRes$n["n"]}}
  \item{\dots}{
}
}



\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{kinrespGrowthphaseExperiment}}
,\code{\link{getUnlimitedGrowthData.kinresp}}
,\code{\link{twKinresp}}}

