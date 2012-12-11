\name{coefKinrespNormStart}
\alias{coefKinrespNormStart}
\title{coefKinrespNormStart}
\description{Transform microbial coefficients to normal scale.}
\usage{coefKinrespNormStart(tmp.coef)}
\arguments{
  \item{tmp.coef}{starting parameters of kinetic respiration, supplied to \code{\link{coefKinresp.default}} (maybe beta)}
}
\details{Replaces \itemize{
\item{x0 by x0l = log(x0)} 
\item{mumax by mumaxl=log(mumax)} and 
\item{r0 by r0l=logit(r0l)}
}
Values at normalized scale are used as starting values for model fits.}
\value{named vector of coefficients (x0l,r0l,mumaxl) }

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{coefKinresp.default}}
,\code{\link{twKinresp}}}

