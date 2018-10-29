\name{mindist}
\alias{mindist}
\title{Mindist measure}
\description{Compute the \code{mindist} criterion (also called maximin)}
\usage{mindist(design)}

\arguments{
  \item{design}{a matrix (or a data.frame) representing the design of experiments in the unit cube [0,1]\eqn{^d}. If this last condition is not fulfilled, a transformation into [0,1]\eqn{^{d}} is applied before the computation of the criteria.}
}

\details{The mindist criterion is defined by
\deqn{mindist= \min_{x_{i}\in X} \left( \gamma_{i} \right)}{mindist = min (g_1, ... g_n)}
where \eqn{\gamma_{i}}{g_i} is the minimal distance between the point \eqn{x_{i}}{x_i}
and the other points \eqn{x_{k}}{x_k} of the \code{design}.

A higher value corresponds to a more regular scaterring of design points.
}

\value{A real number equal to the value of the mindist criterion for the \code{design}.}

\references{

Gunzburer M., Burkdart J. (2004), Uniformity measures for point samples in hypercubes, \url{https://people.sc.fsu.edu/~jburkardt/}.

Jonshon M.E., Moore L.M. and Ylvisaker D. (1990), Minmax and maximin distance designs, \emph{J. of Statis. Planning and Inference}, 26, 131-148.

Chen V.C.P., Tsui K.L., Barton R.R. and Allen J.K. (2003), A review of design and modeling in computer experiments, \emph{Handbook of Statistics}, 22, 231-261.
}

\author{J. Franco}

\seealso{ 
other distance criteria like \code{\link{meshRatio}} and \code{\link{phiP}}, 
discrepancy measures provided by \code{\link{discrepancyCriteria}}.}

\examples{
dimension <- 2
n <- 40
X <- matrix(runif(n*dimension), n, dimension)
mindist(X)
}

\keyword{design}

