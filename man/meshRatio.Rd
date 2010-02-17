\name{meshRatio}
\alias{meshRatio}
\title{MeshRatio measure}
\description{ The \code{meshRatio} criterion is the ratio between the
maximimum and the minimum distance between two points of the experimental design.
}
\usage{meshRatio(design)}

\arguments{
  \item{design}{a matrix (or a data.frame) representing the design of experiments in the unit cube [0,1]\eqn{^d}. If this last condition is not fulfilled, a transformation into [0,1]\eqn{^{d}} is applied before the computation of the criteria.}
}

\details{The \code{meshRatio} criterion is defined by
\deqn{meshRatio=\frac{\max_{1\leq i \leq n} \gamma_{i}}{\min_{1\leq i \leq n}
\gamma_{i}}}{meshRatio(\code{design})=max(g_1, ..., g_n)/min(g_1, ..., g_n)}
where \eqn{\gamma_{i}}{g_i} denotes the minimal distance between the point
\eqn{x_{i}}{x_i} and the other points of the \code{design}.

Note that for a regular mesh, \code{meshRatio}=1.
}
\value{A real number equal to the value of the meshRatio criterion for the \code{design}.}

\references{Gunzburer M. and Burkdart J. (2004) \emph{Uniformity measures for
point samples in hypercubes}
\url{http://people.sc.fsu.edu/~burkardt/pdf/ptmeas.pdf}.}

\author{J. Franco}

\seealso{ other distance criteria like \code{\link{meshRatio}} and \code{\link{mindist}}.

discrepancy measures provided by \code{\link{discrepancyCriteria}}.}

\examples{
dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
meshRatio(X)
}

\keyword{ design }
