\name{CHE1950f}
\alias{CHE1950f}
\docType{data}
\title{
Population and Mortality Data for Swiss Females in 1950
}
\description{
Age-specific deaths, exposures and rates for the Swiss female population in 1950 from the Human Mortality Database.
}
\usage{data("CHE1950f")}
\format{
  A data frame containing the following variables.
  \describe{
    \item{\code{x}}{vector of ages.}
    \item{\code{d}}{vector of death counts.}
    \item{\code{n}}{vector of population exposures.}
    \item{\code{m}}{vector of death rates.}
  }
}
\details{
Data retrieved with the \code{\link{HMD2MH}} function and censored at age 90.}
\source{
Human Mortality Database \url{www.mortality.org}}
\references{
Human Mortality Database. University of California, Berkeley (USA), and Max Planck Institute for Demographic Research (Germany). Available at \url{www.mortality.org} or \url{www.humanmortality.de} (data downloaded in April 2016).}
\examples{
data(CHE1950f)
plot(CHE1950f$x, log(CHE1950f$m), pch = 16, col = 8)}
\seealso{
\code{\link{HMD2MH}}
}
\keyword{datasets}
