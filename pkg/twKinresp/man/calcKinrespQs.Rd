\name{calcKinrespQs}
\alias{calcKinrespQs}
\title{calcKinrespQs}
\description{Calculate respiration components.}
\usage{calcKinrespQs(mumax, YCO2 = 1.5, lambda = 0.9)}
\arguments{
  \item{mumax}{
}
  \item{YCO2}{Yield per respired CO2}
  \item{lambda}{
}
}

\value{list with components \describe{
\item{Q}{maximum substrate uptake for coupled respiration}
\item{Qu}{maximum uncoupled specific respiration}
\item{QT}{Total specific respiration Q_T}
\item{Y}{Yield: biomass/consumed substrate}
}}

\author{Thomas Wutzler <thomas.wutzler@web.de>}



\seealso{\code{\link{twKinresp}}}

