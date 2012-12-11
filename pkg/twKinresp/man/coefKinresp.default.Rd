\name{coefKinresp.default}
\alias{coefKinresp.default}
\title{coefKinresp default}
\description{Check the microbial coefficients and translate to original microbial scale.}
\usage{\method{coefKinresp}{default}(tmp.coef, ...)}
\arguments{
  \item{tmp.coef}{coefficients of model fit \code{coef(model)} (see details)}
  \item{\dots}{
}
}
\details{\describe{\item{Functions for acccessing microbial parameters, and their uncertainty bounds.}{
\itemize{
\item{ Mode, Median, Mean and confidence bounds of microbial fits: \code{\link{kinrespParDist.gnls}}}
\item{ Microbial parameters of single fit: this method (obtained by \code{coef(\link{kinrespGrowthphaseReplicate})} or \code{fixef(\link{fitKinrespExperiment})}) }
\item{ Microbial parameters of replicate fits: \code{\link{coefKinrespMatrix}} (obtained by \code{coef(\link{fitKinrespExperiment})}) }
\item{ Microbial parameters of a list of fits: \code{\link{coefKinresp.kinrespList}} (obtained by \code{\link{kinrespGrowthphaseExperiment}})}
\item{ 95\% confidenc interval of microbial parameters of a single fit: \code{\link{confintKinresp}} }
\item{ Microbial parameters from beta-form of model fit: \code{\link{calcKinrespCoef}} }
}
}}	

\describe{\item{Functions for translating between normal and original scale.}{
\itemize{
\item{ normalized to original scale: all the methods above. }
\item{ original scale to normalized scale: \code{\link{coefKinrespNormStart}} }
}
}}	

\describe{\item{ \code{coef(model)} }{
Models are fitted in various forms differing by used coefficients. 
This method regognizes and translates coefficients of the following forms \itemize{
\item{ beta-form at original scale and at normalized scale}
\item{ beta-form excluding beta0 (fixed to 0 see \code{\link{calcKinrespCoef}}) }
\item{ microbial form fit at original and transformed scale.}
\item{ microbial form fit with excluding r0 (assuming r0=1)}
}
}}}


\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{twKinresp}}}
\examples{FALSE}
