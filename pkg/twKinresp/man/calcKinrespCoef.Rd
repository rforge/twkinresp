\name{calcKinrespCoef}
\alias{calcKinrespCoef}
\title{calcKinrespCoef}
\description{Calculating microbial parameters from beta fit.}
\usage{calcKinrespCoef(tmp.coef, lambda = 0.9, YCO2 = 1.5, 
    cf95 = NA)}
\arguments{
  \item{tmp.coef}{named numeric vector "beta0" \dots "beta2"}
  \item{lambda}{Ratio of growth associated (coupled) specific}
  \item{YCO2}{Ratio of assimilated carbon per respired carbon.}
  \item{cf95}{Confidence intervals of the transformed beta coefficients.}
}



\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{fitKinrespBetaReplicate}}}

