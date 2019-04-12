\name{modelKinrespMic}
\alias{modelKinrespMic}
\title{modelKinrespMic}
\description{Calculate respiration for time, using microbial parameters}
\usage{modelKinrespMic(time, param, lambda = 0.9, 
    YCO2 = 1.5)}
\arguments{
  \item{time}{
}
  \item{param}{named numeric vector ("x0","r0","mumax")}
  \item{lambda}{Ratio of growth associated (coupled) specific}
  \item{YCO2}{Ratio of assimilated carbon per respired carbon.}
}

\value{\code{x0*(1-r0)*(1/lambda-1)*mumax/YCO2 + x0*r0*1/lambda*mumax/YCO2 *
exp(mumax*time)}}

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{twKinresp}}}

