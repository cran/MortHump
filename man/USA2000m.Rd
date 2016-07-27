\name{USA2000m}
\alias{USA2000m}
\docType{data}
\title{
Population and Mortality Data by cause of death for American Males in 2000
}
\description{
Age- and cause-specific death counts and rates for 102 causes of death, population exposures, age labels, age intervals and mid-points, and cause-of-death labels, for the American male population in 2000.}
\usage{data("USA2000m")}
\format{
A list of seven elements containing age- and cause-specific death counts (\code{dxc}) and rates (\code{mxc}), population exposures (\code{nx}), age labels (\code{age}), age intervals (\code{inter}) and mid-points (\code{x}), and cause-of-death labels (\code{lab}), for the American male population in 2000. Data retrieved from the Human Cause-of-Death Database. The mxc and the dxc data frames are structured with ages (five-year age groups) in rows and causes (102 causes) in columns.
}
\details{
Data generated with the \code{\link{HCD2MH}} function.
}
\source{
Human Cause-of-Death Database \url{http://www.causesofdeath.org}
}
\references{
Human Cause-of-Death Database. French Institute for Demographic Studies (France) and Max Planck Institute for Demographic Research (Germany). Available at \url{www.causeofdeath.org} (data downloaded in May 2016).
}
\examples{
data(USA2000m)
matplot(x = USA2000m$x, y = log(USA2000m$mxc[,-1]), type = "l", lty = 1, las = 1)
}
\keyword{datasets}
