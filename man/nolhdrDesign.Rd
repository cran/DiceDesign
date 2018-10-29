\name{nolhdrDesign}
\alias{nolhdrDesign}

\title{De Rainville's Nearly Orthogonal Latin Hypercube Designs}
\description{
     This function generates a NOLH design of dimension 2 to 29 and normalizes it to the selected range. From 2 to 7 input variables, the design is extracted from Cioppa's \code{\link{NOLHdesigns}} list and from 8 to 29 input variables it is extracted from De Rainville's \code{\link{NOLHDRdesigns}} list.
      
}
\usage{
nolhdrDesign(dimension, range = c(0, 1))
}

\arguments{
  \item{dimension }{ number of input variables }
  \item{range     }{ the scale (min and max) of the inputs. Range (0, 0) and (1, 1) are special cases and call integer ranges \eqn{(-m, m)} and \eqn{(0, 2m)}. See the examples}
}

\value{
A list with components:
  \item{n              	}{the number of lines/experiments}
  \item{dimension      	}{the number of columns/input variables}
  \item{design  		}{the design of experiments}
}

\author{T.M. Cioppa and F.-M. De Rainville for the designs. P. Kiener for the R code.}

\seealso{
De Rainville's list \code{\link{NOLHDRdesigns}}. Other NOLH or OLH designs: \code{\link{nolhDesign}}, \code{\link{olhDesign}}.
}

\examples{
## Classical normalizations
nolhdrDesign(8, range = c(1, 1))
nolhdrDesign(8, range = c(0, 1))
nolhdrDesign(8, range = c(0, 0))
nolhdrDesign(8, range = c(-1, 1))

## Change the dimnames, adjust to range (-10, 10) and round to 2 digits
xDRDN(nolhdrDesign(8), letter = "T", dgts = 2, range = c(-10, 10))

## A list of designs
lapply(5:9, function(n) nolhdrDesign(n, range = c(-1, 1))$design)
}

\keyword{design}

