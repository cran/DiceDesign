\name{NOLHDRdesigns}
\alias{NOLHDRdesigns}
\docType{data}

\title{List of De Rainville's Nearly Orthogonal Latin Hypercubes designs}

\description{  
  A list of the NOLH designs for 8 to 29 input variables proposed by De Rainville in 2012. These designs are said to be an improvment of Cioppa's NOLH designs as they have the same structure but better dispersion measures like the discrepancy. 
  
  This list combines the csv files published by De Rainville (see Source), centered and normalized to integer values. It is used internally by the function \code{\link{nolhdrDesign}} which provides various normalizations.

}

\usage{NOLHDRdesigns}

\format{
  A list of 22 matrices representing designs of experiments for 8 to 29 input variables: 
  \describe{
    \item{\code{nolhdr08 to nolhdr11: }}{8, 9, 10, 11 input variables, 33 experiments.}   
    \item{\code{nolhdr12 to nolhdr16: }}{12, 13, 14, 15, 16 input variables, 65 experiments.}
    \item{\code{nolhdr17 to nolhdr22: }}{17, 18, 19, 20, 21, 22 input variables, 129 experiments.}
    \item{\code{nolhdr23 to nolhdr29: }}{23, 24, 25, 26, 27, 28, 29 input variables, 257 experiments.}
  }
}

\author{F.-M. De Rainville for the designs. P. Kiener for the R code.}

\source{ 
Main website: \url{http://qrand.gel.ulaval.ca/}

The python source code: \url{https://github.com/fmder/pynolh/} 

The python package: \url{https://pypi.org/project/pynolh/}.
}

\references{
De Rainville F.-M., Gagne C., Teytaud O., Laurendeau D. (2012). Evolutionary optimization of low-discrepancy sequences. ACM Transactions on Modeling and Computer Simulation (TOMACS), 22(2), 9. 

Cioppa T.M., Lucas T.W. (2007). Efficient nearly orthogonal and space-filling Latin hypercubes. Technometrics 49, 45-55.
} 

\seealso{
The main function \code{\link{nolhdrDesign}}. Cioppa's NOLH design list: \code{\link{NOLHdesigns}}.
}

\examples{
## data(NOLHDRdesigns)

## all matrices 
names(NOLHDRdesigns)
lapply(NOLHDRdesigns, tail, 2)

## The first matrix/design
NOLHDRdesigns[["nolhdr08"]]

}

\keyword{datasets}

