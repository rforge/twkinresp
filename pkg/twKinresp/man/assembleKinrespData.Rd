\name{assembleKinrespData}
\alias{assembleKinrespData}
\title{assembleKinrespData}
\description{Assemble data into form used by the kinresp package.}
\usage{
assembleKinrespData(time, resp, replicate = 1, experiment = "", suite = "")
}
\arguments{
  \item{time}{Time of respriation measurement}
  \item{resp}{Measured respiration}
  \item{replicate}{Identifier of the replicate}
  \item{experiment}{Identifier of the experiment: one replicated time series}
  \item{suite}{Identifier of the suite: Set of experimental condtions}
}
\details{All the parameters need to be vectors of the same length.
The columens experiment and suite are provided to use mixed models to better constrain
those parameters that do not change between experiments.}
\value{data.frame with colums corresponding to the parameters}

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{setFactorKinrespData}}}

