\name{unif.test.statistic}
\alias{unif.test.statistic}
\title{Statistic of some uniformity tests}
\description{Computes the statistic of a uniformity test (see available tests below).}
	
\usage{
unif.test.statistic(x, type, transform=NULL)
}

\arguments{
	\item{x}{a vector containing the sample values.}
	\item{type}{a character indicating which test is used. The choices are the following: "greenwood", "qm" (for Quesenberry-Miller), "ks" (Kolmogorov-Smirnov), "cvm" (Cramer-Von Mises) and "V" (D+ + D- from Kolmogorov-Smirnov).}
		\item{transform}{ an optional character indicating what type of transformation should be applied before testing uniformity. Default is NULL.}

}

\value{A real number equal to the statistic of the specified test.}

\references{ 
	D Agostino R.B., Stephens M.A. (1986), Goodness-of-fit techniques, CRC Press, New York.
}

\author{O. Roustant}

\seealso{\code{\link{unif.test.quantile}}, \code{\link{rss2d}}}

\keyword{design}
