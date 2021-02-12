\name{faureprimeDesign}
\alias{faureprimeDesign}
\title{A special case of the low discrepancy Faure sequence}
\description{
Generate a Faure sequence with \eqn{n=p^u-1}{n=p^u-1} experiments in [0,1]\eqn{^d}{^d} 
or other domains (see the details) where \eqn{p}{p} is the first prime number equal 
or larger than \eqn{d}{d} and \eqn{u}{u} is an exponent, usually 2.
}
\usage{
faureprimeDesign(dimension, u = 2, range = c(0, -1))
}
\arguments{
  \item{dimension}{the number of variables (< 199)}
  \item{u        }{the exponent applied to the prime number}
  \item{range    }{the scale (min and max) of the inputs. See the details for the six predefined ranges.}
}
\details{
This is a special case of \code{\link{runif.faure}} where the number of 
generated points depends exclusively on the dimension and the selected exponent. 
For the exponent \eqn{u=2}{u=2}, the design is orthogonal and has resolution
4. It is a perfect grid \eqn{(p-1)(p+1)}{(p-1)(p+1)} on each pair of 
variables where \eqn{p}{p} is the first prime number equal or 
larger than the dimension \eqn{d}{d}. 


Six domain ranges are predefined and cover most applications:
\itemize{
  \item{\code{c(0, 0)}  corresponds to \eqn{[0, n]^d}{[0, n]^d}.}
  \item{\code{c(1, 1)}  corresponds to \eqn{[1-n, n-1]^d = [2-p^u, p^u -2]^d}{[1-n, n-1]^d = [2-p^u, p^u -2]^d}.}
  \item{\code{c(0, 1)}  corresponds to \eqn{[0, 1]^d}{[0, 1]^d}.}
  \item{\code{c(0,-1)}  corresponds to \eqn{[p^{-u}, 1-p^{-u}]^d}{[p^{-u}, 1-p^{-u}]^d}.}
  \item{\code{c(-1,-1)} corresponds to \eqn{[-1+2p^{-u}, 1-2p^{-u}]^d}{[-1+2p^{-u}, 1-2p^{-u}]^d}.}
  \item{\code{c(-1, 1)} corresponds to \eqn{[-1, 1]^d}{[-1, 1]^d}.}
}
}

\value{
\code{faureprimeDesign} returns a list with the following components:
\itemize{
  \item{design: the design of experiments}
  \item{n: the number of experiments}
  \item{dimension: the dimension}
  \item{prime: the prime number}
  \item{u: the exponent}
}
}

\author{P. Kiener}

\references{
Faure H. (1982), Discrepance de suites associees a un systeme de numeration 
(en dimension s), \emph{Acta Arith.}, 41, 337-351.

Owen A.B. (2020), On dropping the first Sobol point, 
\url{https://arxiv.org/abs/2008.08051}.
}

\examples{
## Range c(0,-1) returns the design produced by runif.faure()
plan1 <- runif.faure(n = 24, dimension = 5)$design ; plan1
plan2 <- faureprimeDesign(dimension = 5, range = c(0,-1))$design ; plan2
all.equal(plan1, plan2, tolerance = 1e-15)

## Range c(0,0) returns the original sequence of integers.
## The first (p-1) lines are on the first diagonal.
## The remaining lines are LHSs grouped in p-1 blocks of p rows.
d <- p <- 5
plan <- faureprimeDesign(dimension = d, range = c(0,0))$design ; plan
apply(plan, 2, sort)

## A regular grid (p-1)x(p+1) rotated by a small angle
pairs(plan)

plot(plan[,1], plan[,2], las = 1)
points(plan[1:(p-1),1], plan[1:(p-1),2], pch = 17, cex = 1.6)
abline(v = plan[1:(p-1),1], col = 4)

## Designs of dimensions 24x5 in various ranges
lstrg <- list(p0p0 = c(0,0),  p1p1 = c(1,1),   p0p1 = c(0,1),
              p0m1 = c(0,-1), m1m1 = c(-1,-1), m1p1 = c(-1,1))
lst <- lapply(lstrg, function(rg) faureprimeDesign(
              dimension = 5, u = 2, range = rg)$design)
lapply(lst, tail)
sapply(lst, range)

## The odd designs (p1m1, m1m1, m1p1) are orthogonal and have resolution 4.
library(lattice)

mat <- lst$m1m1 ; colnames(mat) <- LETTERS[1:5]
fml <- ~ (A+B+C+D+E)^2+I(A^2)+I(B^2)+I(C^2)+I(D^2)+I(E^2)
mmm <- model.matrix(fml, data = as.data.frame(mat))[,-1] ; tail(mmm)
cmm <- round(cov2cor(crossprod(mmm)), 3) ; cmm
lattice::levelplot(cmm[, ncol(cmm):1], at = seq(-1, 1, length.out = 10),
                   col.regions = rev(grDevices::hcl.colors(9, "PuOr")))
}

\keyword{design}

