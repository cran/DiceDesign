\name{rss2d}
\alias{rss2d}
\title{2D graphical tool for defect detection of Space-Filling Designs.}
\description{
For a 2-dimensional design, the 2D radial scanning statistic (RSS) scans angularly the domain. In each direction, it compares the distribution of projected points to their theoretical distribution under the assumption that all design points are drawn from uniform distribution. For a d-dimensional design, all pairs of dimensions are scanned.The RSS detects the defects of low discrepancy sequences or orthogonal arrays, and can be used for selecting space-filling designs.}
	
\usage{
rss2d(design, lower, upper, gof.test.type="greenwood", gof.test.stat=NULL, transform=NULL, n.angle=360, graphics=1, trace=TRUE, lines.lwd = 1, lines.lty = "dotted", ...)
}

\arguments{
	\item{design}{a matrix or data.frame containing the d-dimensional design of experiments. The row no. i contains the values of the d input variables corresponding to simulation no. i}
	\item{lower}{the domain lower boundaries.}
	\item{upper}{the domain upper boundaries.}
	\item{gof.test.type}{ an optional character indicating the kind of statistical test to be used to test the goodness-of-fit of the design projections to their theoretical distribution. Several tests are available, see \code{\link{unif.test.statistic}}. Default is "greenwood".}
	\item{gof.test.stat}{ an optional number equal to the goodness-of-fit statistic at level \code{5\%}. Default is the modified test statistic for fully specified distribution (see details below).}
	\item{transform}{ an optional character indicating what type of transformation should be applied before testing uniformity. Only one choice available "spacings", that lead to over-detection. Default - and recommended - is \code{NULL}.}
	\item{n.angle}{ an optional number indicating the number of angles used. Default is \code{360} corresponding to a 0.5-degree discretization step. Note that the RSS curve is continuous.}
	\item{graphics}{ an optional integer indicating whether a graph should be produced. If negative, no graph is produced. If superior to 2, the RSS curve only is plotted in the worst 2D coordinate subspace (corr. to the worst value of statistic). If 1 (default), the design is also added, with its projections onto the worst (oblique) axis.}
	\item{trace}{ an optional boolean. Turn it to \code{FALSE} if you want no verbosity.}
	\item{lines.lwd}{optional number specifying the width of the straight lines involved in the graphical outputs (axis and projections)}
	\item{lines.lty}{optional character string specifying the type of the straight lines involved in the graphical outputs (axis and projections)}
	\item{...}{ optional graphical parameters of plot function, to draw the RSS curve.}
}


\value{a list with components:
	\item{global.stat}{a matrix containing the values of the global statistic (equal to the maximum of statistic values over the RSS curve) for all pairs of dimensions.}
	\item{worst.case}{the worst pair of dimensions, that is the one that gives the worst value of \code{global.stat}.}
	\item{worst.dir}{the worst direction, that is the one that gives the worst value of the global statistic in the coordinate plane defined by \code{worst.case}.}
	\item{stat}{a vector of length \code{n.angle} containing the statitic values for each angle, in the coordinate plane defined by \code{worst.case}.}
	\item{angle}{a vector of length \code{n.angle} containing the corresponding angles used.}
	\item{curve}{a \code{(2*n.angle)x2} matrix containing the discretized RSS curve.}
	\item{gof.test.stat}{the threshold at significance level 0.05 for the specified goodness-of-fit statistic. It is equal to the radius of the circle superimposed on the RSS figure.}
}

\references{ 
	Roustant O., Franco J., Carraro L., Jourdan A. (2010), A radial scanning statistic for selecting space-filling designs in computer experiments, MODA-9 proceedings, www.emse.fr/~roustant
	
	D Agostino R.B., Stephens M.A. (1986), Goodness-of-fit techniques, CRC Press, New York.
}

\author{O. Roustant}

\seealso{ \code{\link{unif.test.statistic}}, \code{\link{unif.test.quantile}}, \code{\link{rss3d}} }

\keyword{design}

\examples{
# Detection of defects of Sobol designs

# requires randtoolbox package
library(randtoolbox)

	# in 2D
rss <- rss2d(design=sobol(n=20, dim=2), lower=c(0,0), upper=c(1,1), type="l", col="red")

	# in 8D. All pairs of dimensions are tried to detect the worst defect (according to the specified goodness-of-fit statistic).
d <- 8
n <- 10*d
rss <- rss2d(design=sobol(n=n, dim=d), lower=rep(0,d), upper=rep(1,d), type="l", col="red")

# avoid this defect with scrambling ?
#    1. Faure-Tezuka scrambling (type "?sobol" for more details and options)
rss <- rss2d(design=sobol(n=n, dim=d, scrambling=2), lower=rep(0,d), upper=rep(1,d), type="l", col="red")
#    2. Owen scrambling
rss <- rss2d(design=sobol(n=n, dim=d, scrambling=1), lower=rep(0,d), upper=rep(1,d), type="l", col="red")
}

