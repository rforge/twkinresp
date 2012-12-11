\name{kinrespModel1}
\alias{kinrespModel1}
\title{kinrespModel1}
\description{Respiration for given microbial parameters at given time. }
\usage{kinrespModel1(x0, r0, mumax, time, lambda = 0.9, YCO2 = 1.5)}
\arguments{
  \item{x0}{initial microbial biomass (numeric scalar)}
  \item{r0}{initial microbial activity (numeric scalar)}
  \item{mumax}{maximum growth rate (numeric scalar)}
  \item{time}{time (numeric vector)}
  \item{lambda}{Ratio of growth associated (coupled) specific respiration to total specific respiration. Usually 0.9.}
  \item{YCO2}{Ratio of assimilated carbon per respired carbon.  Usually 1.5.}
}

\value{respiration at given time points (numeric vector)}

\author{Thomas Wutzler <thomas.wutzler@web.de>}




\examples{
data(respWutzler10)
ds <- subset(respWutzler10, replicate==1 & experiment ==  1 & suite == 'Face' & time < 25)
resp0 <- kinrespModel1( x0=118, r0=0.029, mumax=0.21, time=ds$time)
plot( resp~time, data=ds)
lines(resp0~ds$time)
}
