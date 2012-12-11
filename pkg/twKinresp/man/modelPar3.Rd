\name{modelPar3}
\alias{modelPar3}
\title{modelPar3}
\description{Translating parameter vector to arguments of Model.}
\usage{modelPar3(parOpt, fModel, parDefault = parOpt, ...)}
\arguments{
  \item{parOpt}{parameters: named numeric vector (x0,r0,mumax)}
  \item{fModel}{model of the form \code{function(x0,r0,mumax,...)}}
  \item{parDefault}{if parOpt contains only a subset of (x0,r0,mumax),
e.g. when optimizing only a subset
parDefault can be used to specify the other values}
  \item{\dots}{further arguments to fModel}
}

\value{result of fModel}

\author{Thomas Wutzler <thomas.wutzler@web.de>}





