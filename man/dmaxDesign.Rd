\name{dmaxDesign}
\alias{dmaxDesign}
\title{ Maximum Entropy Designs }
\description{
  Space-Filling Designs with \eqn{n}{n} experiments based on covariance matrix in
  [0,1]\eqn{^d}{^d}.
}
\usage{
dmaxDesign(n, dimension, range, niter_max=1000, seed=NULL)
}

\arguments{
  \item{n        }{ number of experiments }
  \item{dimension}{ number of variables }
  \item{range    }{ range of variogram }
  \item{niter_max}{ number of iterations}
    \item{seed	 }{ seed used to generate uniform design}
}

\details{Maximum entropy design is a kind of optimal design based
on Shannon's definition of entropy as the amount of information.
Originally, maximum entropy sampling was proposed by Shewry and
Wynn (1987). The goal of the design is to maximize the entropy
defined as the determinant of the correlation matrix using a Fedorov-Mitchell
exchange algorithm.

The spatial correlation matrix is defined by \eqn{C=\left( \rho_{ij}
\right)}{C=(r[i,j])}:
\tabular{ll}{
	\eqn{\rho_{ij}=1-\gamma\left(h_{ij}\right)}{r[i,j]=1-gamma(h[i,j])} \tab if \eqn{h_{ij}\leq a}{h[i,j] <=a}, \cr
	\eqn{\rho_{ij}=0}{r[i,j]=0} \tab if \eqn{h_{ij}>a}{if h_{ij}>a},}
where \eqn{h_{ij}}{h[i,j]} is the distance between \eqn{x_{i}}{xi} and
\eqn{x_{j}}{xj}, \eqn{a}{a} denotes the range of the variogram and
\eqn{\gamma}{gamma} is a spherical variogram:
\deqn{\gamma(h) = 1.5 \frac{h}{a}- 0.5\left(\frac{h}{a}\right)^3
\textnormal{ for } h \leq a}{gamma(h)= 1.5*(h/a)- 0.5*(h/a)^3 for h>=a.}
}

\value{A list with components:
  \item{n              	}{the number of points}
  \item{design			}{the design of experiments}
  \item{dimension      	}{the number of variables}
  \item{range          	}{the range of the variogram}
  \item{niter_mx		}{the number of iterations}
  \item{design_init		}{the initial distribution}
  \item{det_init       	}{the value of the determinant for the initial distribution}
  \item{det_end        	}{the value of the determinant at the end of the procedure}
  \item{seed			}{the value of the seed}
}

\references{Currin C., Mitchell T., Morris M. and Ylvisaker D. (1991)
\emph{Bayesian Prediction of Deterministic Functions With Applications
to the Design and Analysis of Computer Experiments}, American
Statistical Association, \bold{86}, 416, 953-963.

Shewry, M. C. and Wynn and H. P. (1987) \emph{Maximum entropy sampling},
Journal of Applied Statistics 14, 165-170.}

\author{J. Franco}

\examples{
n <- 20
dimension <- 2
range <-0.9
niter_max <- 200
out <- dmaxDesign(n, dimension, range, niter_max)

## Change the dimnames, adjust to range (-10, 10) and round to 2 digits
xDRDN(out, letter = "T", dgts = 2, range = c(-10, 10))
}

\keyword{design}

