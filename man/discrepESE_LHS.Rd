\name{discrepESE_LHS}
\alias{discrepESE_LHS}

\title{Enhanced Stochastic Evolutionnary (ESE) algorithm for Latin Hypercube Sample (LHS) optimization via L2-discrepancy criteria}
\description{The objective is to produce low-discrepancy LHS. ESE is a powerful genetic algorithm to produce space-filling designs. It has been adapted here to main discrepancy criteria.}

\usage{
discrepESE_LHS(design, T0=0.005*discrepancyCriteria(design,type='C2')[[1]], 
inner_it=100, J=50, it=2, criterion="C2")
}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments.}
  \item{T0}{The initial temperature  of the ESE algorithm}
  \item{inner_it}{The number of iterations for inner loop}
  \item{J}{The number of new proposed LHS inside the inner loop}
  \item{it}{ The number of iterations for outer loop}
  \item{criterion}{The criterion to be optimized. One can choose three different L2-discrepancies: the C2 (centered) discrepancy ("C2"), the L2-star discrepancy ("L2star") and the W2 (wrap-around) discrepancy ("W2")}
  }

\details{This function implements a stochastic algorithm (ESE) to produce optimized LHS. It is based on Jin et al works (2005). Here, it has been adapted to some discrepancy criteria taking into account new ideas about the revaluations of discrepancy value after a LHS elementary perturbation (in order to avoid computing all terms in the discrepancy formulas).}


\value{A list containing:
  \item{InitialDesign}{the starting design}
  \item{T0}{the initial temperature of the ESE algorithm}
  \item{inner_it}{the number of iterations for inner loop}
  \item{J}{the number of new proposed LHS inside the inner loop}
  \item{it}{the number of iterations for outer loop}
   \item{criterion}{the criterion to be optimized}
  \item{design}{the matrix of the final design (low-discrepancy LHS)}
  \item{critValues}{vector of criterion values along the iterations}
  \item{tempValues}{vector of temperature values along the iterations}
  \item{probaValues}{vector of acceptation probability values along the iterations}
}

\references{

Damblin G., Couplet M., and Iooss B. (2013). Numerical studies of space filling designs: optimization of Latin Hypercube Samples and subprojection properties, \emph{Journal of Simulation,} 7:276-289, 2013.
\url{http://www.gdr-mascotnum.fr/doku.php?id=iooss1}

M. Morris and J. Mitchell (1995) Exploratory designs for computational experiments. Journal of Statistical Planning and Inference, 43:381-402.

R. Jin, W. Chen and A. Sudjianto (2005) An efficient algorithm for constructing optimal design of computer experiments. Journal of Statistical Planning and Inference, 134:268-287.

}

\author{G.Damblin & B. Iooss}

\seealso{
Latin Hypercube Sample(\code{\link{lhsDesign}}), 
discrepancy criteria(\code{\link{discrepancyCriteria}}), 
geometric criterion (\code{\link{mindist}}\code{\link{phiP}}), 
optimization (\code{\link{maximinSA_LHS}}, \code{\link{maximinESE_LHS}}, \code{\link{discrepSA_LHS}})
}


\examples{
\dontrun{
dimension <- 2
n <- 10
X <- lhsDesign(n, dimension)$design
Xopt <- discrepESE_LHS(X, T0=0.005*discrepancyCriteria(X, type='C2')[[1]], 
                       inner_it=100, J=50, it=2)
plot(Xopt$design)
plot(Xopt$critValues, type="l")
}
}

\keyword{design}

