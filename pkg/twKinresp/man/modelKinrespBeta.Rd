\name{modelKinrespBeta}
\alias{modelKinrespBeta}
\title{modelKinrespBeta}
\description{Calculate respiration for time x based on coefficients beta_i.}
\usage{modelKinrespBeta(x, betai)}
\arguments{
  \item{x}{time}
  \item{betai}{named numeric vector (beta0,beta1,beta2)}
}

\value{\code{beta0 + beta1*exp(beta2*x)}}

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{fitKinrespBetaReplicate}}}

