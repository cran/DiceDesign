\name{phiP}
\alias{phiP}

\title{phiP criterion}
\description{Compute the \eqn{\phi_p} criterion (strongly linked to mindist criterion)}
\usage{phiP(design, p=50)}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments.}
  \item{p}{the "p" in the Lp norm which is taken} 
  }


\details{The \eqn{\phi_p} criterion is defined by the \eqn{L_p} norm of the sum of the inverses of the design inter-point euclidean distances:
\deqn{\phi_{p}=\left[\sum_{i,j=1\ldots N,i<j}\,\,d_{ij}^{-p}\right]^{\frac{1}{p}}}{}

A higher value corresponds to a more regular scaterring of design points.

When \eqn{p} tends to infinity, optimizing a design with \eqn{\phi_p} is equivalent to optimizing a design with \code{mindist}.
}


\value{A real number equal to the value of the \eqn{\phi_p} criterion for the \code{design}.}

\references{

Damblin G., Couplet M., and Iooss B. (2013). Numerical studies of sapce filling designs: optimization of Latin Hypercube Samples and subprojection properties, \emph{Journal of Simulation},  7:276-289, 2013. \url{http://www.gdr-mascotnum.fr/doku.php?id=iooss1}

Fang K.-T., Li R. and Sudjianto A. (2006). Design and Modeling for Computer Experiments, \emph{Chapman & Hall}.

Pronzato, L. and Muller, W. (2012). Design of computer experiments: space filling and beyond, \emph{Statistics and Computing}, 22:681-701.

}

\author{G. Damblin & B.Iooss}

\seealso{geometric criterion (\code{\link{mindist}})} 


\examples{
dimension <- 2
n <- 40
X <- matrix(runif(n*dimension), n, dimension)
phiP(X)
}

\keyword{design}

