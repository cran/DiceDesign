\name{coverage}
\alias{coverage}
\title{Coverage}
\description{Compute the coverage measure}
\usage{coverage(design)}

\arguments{
  \item{design}{a matrix (or a data.frame) representing the design of experiments representing the design of experiments in the unit cube [0,1]\eqn{^d}. If this last condition is not fulfilled, a transformation into [0,1]\eqn{^{d}} is applied before the computation of the criteria.}
}

\details{The coverage criterion is defined by
\deqn{coverage=\frac{1}{\bar{\gamma}} \left[ \frac{1}{n} \sum_{i=1}^{n}
	\left( \gamma_{i} - \bar{\gamma} \right)^2 \right]^{1/2}}{coverage)
	=1/gMean *[ 1/n * [( g_1 - gMean )^2 + ... + (g_n - gMean)^2] ]^(1/2)}
where \eqn{\gamma_{i}}{g_i} is the minimal distance between the point \eqn{x_{i}}{x_i}
and the other points of the \code{design} and \eqn{\bar{\gamma}}{gMean} is
the mean of the \eqn{\gamma_{i}}{g_i}. 

Note that for a regular mesh, \code{cov}=0. Then, a small value of \code{cov} means that the design is close to a regular grid.}
\value{A real number equal to the value of the coverage criterion for the \code{design}.}

\references{Gunzburer M., Burkdart J. (2004) \emph{Uniformity measures for point samples in hypercubes}, \url{https://people.sc.fsu.edu/~jburkardt/}.}

\author{J. Franco}

\seealso{other distance criteria like \code{\link{meshRatio}}, \code{\link{phiP}} and \code{\link{mindist}}.

discrepancy measures provided by \code{\link{discrepancyCriteria}}.}

\examples{
dimension <- 2
n <- 40
X <- matrix(runif(n*dimension), n, dimension)
coverage(X)
}

\keyword{design}

