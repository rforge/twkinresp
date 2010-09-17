\name{fitKinrespExpAll}
\alias{fitKinrespExpAll}
\title{fitKinrespExpAll}
\description{Invoking fitKinrespExp for each experiment within the suite.}
\usage{
fitKinrespExpAll(rd.e, repFits, lambda = 0.9, YCO2 = 1.5, weights = NULL)
}
\arguments{
  \item{rd.e}{dataset with columns suite, experiment, replicate, resp and time, containing only unlimited growth phase.}
  \item{repFits}{Initial coefficients mumax, x0, and r0 for each replicate}
  \item{lambda}{
}
  \item{YCO2}{
}
  \item{weights}{
}
}



\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{twKinresp}}}

