\name{fitKinrespBetaReplicate}
\alias{fitKinrespBetaReplicate}
\title{fitKinrespBetaReplicate}
\description{Fitting the beta-model with estimating the log of beta0, so that exp will be positive.}
\usage{
fitKinrespBetaReplicate(x, y, weights = NULL)
}
\arguments{
  \item{x}{
}
  \item{y}{
}
  \item{weights}{
}
}
\details{\describe{\item{Functions related to fitting the beta-form of the model.}{
\describe{
\item{\code{\link{modelKinrespBeta}} }{ Model form with simple coefficients beta0..beta2.  }
\item{\code{\link{calcKinrespCoef}} }{ Calulate Microbial kinetic parameters from given beta coefficients.  }
\item{\code{\link{coefKinrespBeta.default}} }{ Transform coefficients of Beta-form model from transformed scale to original scale.  }
}
}}}


\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{twKinresp}}
, \code{\link{kinrespGrowthphaseReplicate}}}

