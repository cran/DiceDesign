\name{NOLHdesigns}
\alias{NOLHdesigns}
\docType{data}

\title{List of Cioppa's Nearly Orthogonal Latin Hypercubes designs}
\description{
   A list of the NOLH designs for 2 to 29 input variables proposed by Cioppa in 2007. These designs combine a latin structure, orthogonality between the main terms and the interactions (+ squares) and reduced correlations between the interactions (+ squares). 
  
   This list combines the Excel spreadsheets published by Sanchez (see Source). It is used internally by the function \code{\link{nolhDesign}} which provides various normalizations.
}

\usage{NOLHdesigns}
\format{
  A list of 5 matrices representing designs of experiments for 8 to 29 input variables: 
  \describe{
    \item{\code{nolh2_7: }}{2 to 7 input variables, 17 experiments.}
    \item{\code{nolh8_11: }}{8 to 11 input variables, 33 experiments.}
    \item{\code{nolh12_16: }}{12 to 16 input variables, 65 experiments.}
    \item{\code{nolh17_22: }}{17 to 22 input variables, 129 experiments.}
    \item{\code{nolh23_29: }}{23 to 29 input variables, 257 experiments.}
  }
}

\author{T.M. Cioppa for the designs. P. Kiener for the R code.}

\source{
Sanchez, S. M. (2011). NOLHdesigns in Excel file. Available online at \url{https://nps.edu/web/seed/software-downloads/} 
} 

\references{
Cioppa T.M., Lucas T.W. (2007). Efficient nearly orthogonal and space-filling Latin hypercubes. Technometrics 49, 45-55.

Kleijnen, J.P.C., Sanchez S.M., T.W. Lucas and Cioppa T. M.. A user's guide to the brave new world of designing simulation experiments. INFORMS Journal on Computing 17(3): 263-289.

Ye, K. Q. (1998). Orthogonal Latin hypercubes and their application in computer experiments. J. Amer. Statist. Asso. 93, 1430-
1439.
}

\seealso{
The main function \code{\link{nolhDesign}}. De Rainville's NOLH design list: \code{\link{NOLHDRdesigns}}.
}

\examples{
## data(NOLHdesigns)

## all matrices 
names(NOLHdesigns)
lapply(NOLHDRdesigns, tail, 2)

## The first matrix/design
NOLHdesigns[["nolh2_7"]]
}

\keyword{datasets}

