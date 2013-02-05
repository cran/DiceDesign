\name{wspDesign}
\alias{wspDesign}

\title{WSP algorithm}
\description{The WSP (Wooton, Sergent, Phan-Tan-Luu) algorithm is an iterative algorithm based on suppression of some experiments from an initial design in each step. WSP leads to a space filling design.}
\usage{wspDesign(design,dmin)}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments.}
  \item{dmin}{a minimum bound for mindist value of the final design} 
  }


\details{WSP enables to create a design D which is such that mindist(D)>dmin. However, it cannot assess the number of experiments. Similarly to straussDesign function, WSP is a powerful algorithm to construct space filling designs in high dimension} 


\value{A matrix corresponding to the design of experiments.}

\references{

J. Santiago, M. Claeys-Bruno, M.Sergent (2012). Construction of space filling designs using WSP algorithm for high dimensional spaces,
\emph{Chenometrics and Intelligent Laboratory Systems 113:26-31}

}

\author{G.Damblin & B.Iooss}


\examples{
dimension <- 2
n <- 100
X <- matrix(runif(n*dimension),n,dimension)
m=wspDesign(X,0.1)
plot(m)
}

\keyword{ design }