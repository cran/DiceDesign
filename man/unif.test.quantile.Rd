\name{unif.test.quantile}
\alias{unif.test.quantile}
\title{Quantile of some uniformity tests}
\description{Computes the quantile of a uniformity test at a given significance level (see available tests and levels below).}
	
\usage{
unif.test.quantile(type, n, alpha)
}

\arguments{
	\item{type}{a character indicating which test is used. The choices are the following: "greenwood", "qm" (for Quesenberry-Miller), "ks" (Kolmogorov-Smirnov), "cvm" (Cramer-Von Mises) and "V" (D+ + D- from Kolmogorov-Smirnov).}
	\item{n}{an integer equal to the sample size.}
	\item{alpha}{a real number equal to significance level. At present stage, only four values are available: 0.1, 0.05, 0.025 and 0.01.}
}

\value{A real number equal to the quantile of the specified test at significance level \code{alpha} for \code{n} observations.}

\references{ 
	D Agostino R.B., Stephens M.A. (1986), Goodness-of-fit techniques, CRC Press, New York.
}

\details{
	Modified statistics are used. For \code{alpha = 0.05}, the quantile is (see D Agostino and Stephens, 1986, section 4.4.): \code{1.358/(sqrt(n) + 0.12 + 0.11/sqrt(n))} for Kolmogorov-Smirnov and \code{0.461/(1+1/n) + 0.4/n - 0.6/n^2} for Cramer-von Mises. When the design size is \code{< 20}, the corrected value seems to be a good approximation, but the non asymptotical value should be preferred.
}

\author{O. Roustant}

\seealso{\code{\link{unif.test.statistic}}, \code{\link{rss2d}}, \code{\link{rss3d}}}

\keyword{design}
