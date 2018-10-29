\name{lhsDesign}
\alias{lhsDesign}
\title{Latin Hypercube Designs}
\description{
  Simple (random) Latin Hypercube Design (randomized or centered) with \eqn{n}{n} experiments in [0,1]\eqn{^d}{^d}.
}
\usage{
lhsDesign(n, dimension, randomized=TRUE, seed=NULL)
}

\arguments{
  \item{n         }{ number of experiments }
  \item{dimension }{ number of variables }
  \item{randomized}{ TRUE for randomized LHS; FALSE for centered LHS}
  \item{seed      }{ seed used to generate the random permutations and perturbations}
}

\details{This program builds a Latin Hypercube Design (LHD), also called a Latin Hypercube Sample (LHS), on the space [0,1]\eqn{^d}{^d} (with uniform probability measures). LHD aims at ensuring that each variable has its whole range well scanned: the range of each variable is divided into \eqn{n}{n} equally probable stratas. Each stratum of each variable contains only one point of the LHD. Centered LHD is obtained by choosing for each point the center of the corresponding case, while randomized LHD is obtained by adding random perturbations inside each point case.

Once the sample is generated, the uniform sample from a column can be transformed to any distribution by using the quantile functions.
}

\value{A list with components:
  \item{n              	}{the number of points}
  \item{dimension      	}{the number of variables}
  \item{design  	    }{the design of experiments}
  \item{randomized      }{the type of LHD}
  \item{seed			}{the value of the seed}
}

\references{McKay M., Conover W. and Beckman R. (1979)
\emph{A comparison of three methods for selecting values of input variables in the analysis of output from a computer code}, Technometrics, \bold{21}, 2, 239-245.

Stein M. (1987)
\emph{Large sample properties of simulations using Latin hypercube sampling}, Technometrics, \bold{29}, 143-151.
}

\author{B. Iooss}

\seealso{LHD optimization (\code{\link{maximinSA_LHS}},\code{\link{discrepSA_LHS}},\code{\link{maximinESE_LHS}} ,\code{\link{discrepESE_LHS}})}

\examples{
n <- 20
dimension <- 2
out <- lhsDesign(n, dimension)
out$design

## Change the dimnames, adjust to range (-10, 10) and round to 2 digits
xDRDN(out, letter = "T", dgts = 2, range = c(-10, 10))
}

\keyword{design}

