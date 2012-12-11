\name{fitKinrespReplicate}
\alias{fitKinrespReplicate}
\title{fitKinrespReplicate}
\description{Fitting the kinResp microbial model to time series of one replicate.}
\usage{fitKinrespReplicate(rder.e, lambda = 0.9, YCO2 = 1.5, start = NULL, 
    weights = NULL)}
\arguments{
  \item{rder.e}{Respiration dataset containing columns replicate, resp and time, constrained to unlimited growth phase.}
  \item{lambda}{
}
  \item{YCO2}{
}
  \item{start}{
}
  \item{weights}{
}
}
\details{If the microbial explicit form did not fit, then the beta fit is returned.
The beta fit itself might not contain beta0, because it was constrained to 1.}


\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{twKinresp}}, \code{\link{coefKinresp}} }

