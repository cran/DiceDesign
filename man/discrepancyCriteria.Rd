\name{discrepancyCriteria}
\alias{discrepancyCriteria}

\title{Discrepancy measure}
\description{Compute discrepancy criteria. }
\usage{discrepancyCriteria(design,type='all')}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments.
  The discrepancy criteria are computed for a design in the unit cube [0,1]\eqn{^d}.
  If this condition is not satisfied the design is automatically rescaled.}
  \item{type}{type of discrepancies (single value or vector) to be computed:
	\tabular{ll}{
		\code{'all'} \tab all type of discrepancies (default) \cr
		\code{'L2'} \tab star L2-discrepancy \cr
		\code{'C2'} \tab centered L2-discrepancy \cr
		\code{'M2'} \tab modified L2-discrepancy \cr
		\code{'S2'} \tab symmetric L2-discrepancy \cr
		\code{'W2'} \tab wrap-around L2-discrepancy \cr
	}
  }
}

\details{The discrepancy measures how far a given distribution of points deviates
from a perfectly uniform one. The computed values represent the quantity 
\deqn{\left\| \frac{A\left(J,X\right)}{n} - Vol \left( J \right) \right\|_{2}}{|| 
A (J,X)/n - Vol(J))||_{L2}} where \eqn{A\left(J,X\right)}{A (J,X)} is the number 
of points of the design \eqn{X} falling in a part \eqn{J} of the unit cube depending 
on the type of discrepancy.

Centered L2-discrepancy is computed using the analytical expression done by Hickernell
(1998). The user will refer to Pleming and Manteufel (2005) to have more details about 
the wrap around discrepancy.
}

\value{A real number equal to the centered L2-discrepancy of the \code{design}.}

\references{Fang K.T, Li R. and Sudjianto A. (2006) Design and Modeling for
	Computer Experiments, \emph{Chapman & Hall}.

Franco J. (2008) Planification d'experiences numerique en phase exploratoire pour
la simulation des phenomenes complexes, \emph{PhD thesis, Ecole Nationale Superieure
des Mines de Saint Etienne}.

Hickernell F.J. (1998) A generalized discrepancy and quadrature error bound.
	\emph{Mathematics of Computation}, \bold{67}, 299-322.

Pleming J.B. and Manteufel R.D. (2005) \emph{Replicated Latin Hypercube Sampling},
46th Structures, Structural Dynamics & Materials Conference, 16-21 April 2005, Austin
(Texas) -- AIAA 2005-1819.
}

\author{J. Franco & D. Dupuy}

\seealso{distance criteria (\code{\link{coverage}}, \code{\link{meshRatio}} and
\code{\link{mindist}})}

\examples{
dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
discrepancyCriteria(X)
}

\keyword{ design }

