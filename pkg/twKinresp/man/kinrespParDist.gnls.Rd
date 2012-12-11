\name{kinrespParDist.gnls}
\alias{kinrespParDist.gnls}
\title{kinrespParDist gnls}
\description{Expected value and confidence interval of microbial parameters}
\usage{\method{kinrespParDist}{gnls}(model, ...)}
\arguments{
  \item{model}{the result of the \code{\link{fitKinrespReplicate}} or \code{\link{fitKinrespExperiment}$model}}
  \item{\dots}{currently not used}
}
\details{Both the lognormal and the logitnormal distributions are skewed.
Therefore, the mode, the median, and the expected value differ.
At transformed, i.e. normal scale, they coincide. Simple backtransformation
yield the median at original scale.}
\value{numeric matrix with rows corresponding to variables and columns \itemize{
\item{ \code{mean}: expected value}
\item{ \code{sd}: standard deviation}
\item{ \code{mle}: the mode, i.e. the maximum likelihood estimate}
\item{ \code{median}: the median}
\item{ \code{cf025} and \code{cf975}: 2.5\% and 97.5\% confidence bounds}
\item{ \code{mu} and \code{sigma}: parameters at normal, i.e log or logit transformed scale}
}}

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{coefKinresp.default}}
,\code{\link{twKinresp}}}
\examples{FALSE}
