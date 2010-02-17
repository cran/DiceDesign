\name{runif.faure}
\alias{runif.faure}
\title{Low discrepancy sequence : Faure }
\description{
	Generate a Faure sequence with \eqn{n}{n} experiments in [0,1]\eqn{^d}{^d}.}
\usage{
runif.faure(n, dimension)
}

\arguments{
  \item{n        }{the number of experiments}
  \item{dimension}{the number of variables (<100)}
}

\details{A quasirandom or low discrepancy sequence, such as the
Faure, Halton, Hammersley, Niederreiter or Sobol sequences,
is "less random" than a pseudorandom number sequence, but more
useful for such tasks as approximation of integrals in higher
dimensions, and in global optimization. This is because low
discrepancy sequences tend to sample space "more uniformly"
than random numbers.

see \pkg{randtoolbox} or \pkg{fOptions} packages for other low discrepancy sequences.}

\value{\code{runif.halton} returns a list containing all the
input arguments detailed before, plus the following component:
  \item{design   }{the design of experiments}
}
\references{Faure H. (1982) \emph{Discrépance de suites associées
à un système de numération (en dimension s)}, Acta Arith. 41,
337-351 }

\author{J. Franco}

\examples{
f <- runif.faure(20,2)
plot(f$design,xlim=c(0,1),ylim=c(0,1))
}

\keyword{design }