\name{factDesign}
\alias{factDesign}
\title{Full Factorial Designs}
\description{
Create a factorial design with \eqn{n}{n} = pow(levels,dimension) experiments
in [0,1]\eqn{^d}{^d}.
}

\usage{
factDesign(dimension, levels)
}

\arguments{
  \item{dimension}{an integer given the number of input variables }
  \item{levels   }{an integer given the number of levels}
}

\details{It is possible to take a different number of levels for
any factor. In this case, the argument \code{levels} should be a
vector.
}

\value{\code{factDesign} returns a list containing all the
input arguments detailed before, plus the following components:
  \item{n         }{the number of experiments}
  \item{design    }{the design of experiments}
}
\author{G. Pujol and J. Franco}


\examples{
## First example
g1 <- factDesign(2, 7)
plot(g1$design, xlim=c(0,1), ylim=c(0,1))

## Second example
g2 <- factDesign(2, c(2,7))
plot(g2$design, xlim=c(0,1), ylim=c(0,1))

## Change the dimnames, adjust to range (-10, 10) and round to 2 digits
xDRDN(g1, letter = "T", dgts = 2, range = c(-10, 10))
xDRDN(g2, letter = "T", dgts = 2, range = c(-10, 10))
}

\keyword{design}

