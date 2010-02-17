\name{rss3d}
\alias{rss3d}
\title{3D graphical tool for defect detection of Space-Filling Designs.}
\description{
For a 3-dimensional design, the 3D radial scanning statistic (RSS) scans angularly the domain. In each direction, it compares the distribution of projected points to their theoretical distribution under the assumption that all design points are drawn from uniform distribution. For a d-dimensional design, all triplets of dimensions are scanned. The RSS detects the defects of low discrepancy sequences or orthogonal arrays, and can be used for selecting space-filling designs.}

\usage{
rss3d(design, lower, upper, gof.test.type = "greenwood", gof.test.stat = NULL, transform = NULL, n.angle = 60, graphics = 1, trace = TRUE)
}

\arguments{
	\item{design}{a matrix or data.frame containing the d-dimensional design of experiments. The row no. i contains the values of the d input variables corresponding to simulation no. i}
	\item{lower}{the domain lower boundaries.}
	\item{upper}{the domain upper boundaries.}
	\item{gof.test.type}{an optional character indicating the kind of statistical test to be used to test the goodness-of-fit of the design projections to their theoretical distribution. Several tests are available, see \code{\link{unif.test.statistic}}. Default is "greenwood".}
	\item{gof.test.stat}{an optional number equal to the goodness-of-fit statistic at level \code{5\%}. Default is the modified test statistic for fully specified distribution (see details below).}
	\item{transform}{ an optional character indicating what type of transformation should be applied before testing uniformity. Only one choice available "spacings", that lead to over-detection. Default - and recommended - is \code{NULL}.}
	\item{n.angle}{ an optional number indicating the number of angles used. Default is \code{60} corresponding to a 3-degree discretization step. Note that the RSS surface is continuous.}
	\item{graphics}{ an optional integer indicating whether a graph should be produced. If negative, no graph is produced. Otherwise (default), the design is plotted in the worst 3D coordinate subspace (corr. to the worst value of statistic), with its projections onto the worst (oblique) axis.}
	\item{trace}{ an optional boolean. Turn it to \code{FALSE} if you want no verbosity.}
}

\details{
	The RSS surface is continuous. However for computational purposes, a discretization is used. The default discretization step is tunable with \code{n.angle}.
}


\value{a list with components:
	\item{global.stat}{an array containing the values of the global statistic (equal to the maximum of statistic values over the RSS surface) for all triplets of dimensions.}
	\item{print.out}{the same as \code{global.stat}, but with a user-friendly printing.}
	\item{worst.case}{the worst triplet of dimensions, that is the one that gives the worst value of \code{global.stat}.}
	\item{worst.dir}{the worst direction, that is the one that gives the worst value of the statistic in the coordinate 3D subspace defined by \code{worst.case}.}
	\item{stat}{a matrix of size \code{n.angle*n.angle} containing the statistic values for each angles (spherical coordinates).}
	\item{angle}{a matrix of size \code{n.angle*n.angle} containing the corresponding angles used (spherical coordinates).}
	\item{gof.test.stat}{the threshold at significance level 0.05 for the specified goodness-of-fit statistic.}
}



\references{ 
	Roustant O., Franco J., Carraro L., Jourdan A. (2010), A radial scanning statistic for selecting space-filling designs in computer experiments, MODA-9 proceedings, www.emse.fr/~roustant
	
	D Agostino R.B., Stephens M.A. (1986), Goodness-of-fit techniques, CRC Press, New York.
}

\author{O. Roustant}


\seealso{ \code{\link{unif.test.statistic}}, \code{\link{unif.test.quantile}}, \code{\link{rss2d}}}

\keyword{design}

\examples{
	# An orthogonal array in 3D   
data(OA131)

	# centering the design points of this 7-levels design
OA <- (OA131 + 0.5)/7

	# 2D projections onto coordinate axis 
pairs(OA, xlim=c(0,1), ylim=c(0,1))
	# Now let us look at the 3D properties with the 3D RSS (requires the rgl package)
rss <- rss3d(OA, lower=c(0,0,0), upper=c(1,1,1))
	# The worst direction detected is nearly proportional to (2,-1,2)
	# (type "?OA131" for explanations about this linear orthogonal array)
print(rss$worst.dir)

# Now, scramble this design
# X <- (OA131 + matrix(runif(49*3, 49, 3)))/7
# or load the design obtained this way
data(OA131_scrambled)
OA2 <- OA131_scrambled

rss <- rss2d(OA2, lower=c(0,0,0), upper=c(1,1,1))    # no feature is detected by the 2D RSS
rss <- rss3d(OA2, lower=c(0,0,0), upper=c(1,1,1))	# 4 clusters are detected by the 3D RSS

	
	# Defect detection of 8D Sobol sequences
	# All triplets of dimensions are tried to detect the worst defect (according to the specified goodness-of-fit statistic).
	# requires randtoolbox library to generate the Sobol sequence
\dontrun{
library(randtoolbox)
d <- 8
n <- 10*d
rss <- rss3d(design=sobol(n=n, dim=d), lower=rep(0,d), upper=rep(1,d))}
}

