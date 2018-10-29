\name{olhDesign}
\alias{olhDesign}

\title{Nguyen's Orthogonal Latin Hypercube Designs}
\description{
   Generate the Orthogonal Latin Hypercube (OLH) designs proposed by Nguyen in 2008. These OLHs have a latin structure, an orthogonality between the main terms and the interactions (+ squares) and low correlations between the interactions (+ squares). Very larges matrices can be obtained as the number of input variables and hence the number of lines is unconstrained. When the number of input variables is a power of 2, OLHs have \eqn{d} columns and \eqn{n = 2d + 1} lines (experiments). A vertical truncature is applied when the number of input variables is not a power of 2. Various normalizations can be applied. 
      
}
\usage{
olhDesign(dimension, range = c(0, 1))
}

\arguments{
  \item{dimension }{ number of input variables }
  \item{range     }{ the scale (min and max) of the inputs. Ranges (0, 0) and (1, 1) are special cases and call integer ranges \eqn{(-d, d)} and \eqn{(0, 2d)}. See the examples}
}

\value{
A list with components:
  \item{n              	}{the number of lines/experiments}
  \item{dimension      	}{the number of columns/input variables}
  \item{design  		}{the design of experiments}
}

\author{N.K. Nguyen for the algorithm. P. Kiener for the recursive R code.}

\references{
Nguyen N.K. (2008) \emph{A new class of orthogonal Latinhypercubes}, Statistics and Applications, Volume 6, issues 1 and 2, pp.119-123.
}

\seealso{
Cioppa's and De Rainville's NOLH designs: \code{\link{nolhDesign}}, \code{\link{nolhdrDesign}}.
}

\examples{
## Classical normalizations
olhDesign(4, range = c(0, 0))
olhDesign(4, range = c(1, 1))
olhDesign(4, range = c(0, 1))
olhDesign(4, range = c(-1, 1))

## Change the dimnames, adjust to range (-10, 10) and round to 2 digits
xDRDN(olhDesign(4), letter = "T", dgts = 2, range = c(-10, 10))

## A list of designs
lapply(1:5, function(n) olhDesign(n, range = c(-1, 1))$design)
}

\keyword{design}

