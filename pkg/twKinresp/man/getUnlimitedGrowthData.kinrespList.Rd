\name{getUnlimitedGrowthData.kinrespList}
\alias{getUnlimitedGrowthData.kinrespList}
\title{getUnlimitedGrowthData kinrespList}
\description{Extract the dataset with records of unlimited growth phase.}
\usage{\method{getUnlimitedGrowthData}{kinrespList}(kinrespRes, 
    n = integer(0), ...)}
\arguments{
  \item{kinrespRes}{object of class kinrespList from \code{\link{kinrespGrowthphaseExperiment}}.}
  \item{n}{integer vector with names (character suite_experiment_replicate) speciying the number of records in growth phase,
<< if not given (length zero), the defaults from \code{kinrespRes$n["n"]} are used}
  \item{\dots}{
}
}



\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{kinrespGrowthphaseExperiment}}
,\code{\link{getUnlimitedGrowthData.kinresp}}
,\code{\link{twKinresp}}}

