\name{discrepSA_LHS}
\alias{discrepSA_LHS}

\title{Simulated annealing (SA) routine for Latin Hypercube Sample (LHS) optimization via L2-discrepancy criteria}
\description{The objective is to produce low-discrepancy LHS. SA is an efficient algorithm to produce space-filling designs. It has been adapted here to main discrepancy criteria.}

\usage{
discrepSA_LHS(design, T0=10, c=0.95, it=2000, criterion="C2", profile="GEOM", Imax=100)
}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments}
  \item{T0}{The initial temperature}
  \item{c}{A constant parameter regulating how the temperature goes down}
  \item{it}{The number of iterations}
  \item{criterion}{The criterion to be optimized. One can choose three different L2-discrepancies: the C2 (centered) discrepancy ("C2"), the L2-star discrepancy ("L2star") and the W2 (wrap-around) discrepancy ("W2")}
  \item{profile}{The temperature down-profile, purely geometric called "GEOM", geometrical according to the Morris algorithm called "GEOM_MORRIS" or purely linear called "LINEAR"}
  \item{Imax}{A parameter given only if you choose the Morris down-profile. It adjusts the number of iterations without improvement before a new elementary perturbation}
  }


\details{This function implements a classical routine to produce optimized LHS. It is based on the work of Morris and Mitchell (1995). They have proposed a SA version for LHS optimization according to mindist criterion. Here, it has been adapted to some discrepancy criteria taking in account new ideas about the reevaluations of a discrepancy value after a LHS elementary perturbation (in order to avoid computing all terms in the discrepancy formulas).}


\value{A list containing:
  \item{InitialDesign}{the starting design}
  \item{T0}{the initial temperature of the SA algorithm}
  \item{c}{the constant parameter regulating how the temperature goes down}
  \item{it}{the number of iterations}
   \item{criterion}{the criterion to be optimized}
  \item{profile}{the temperature down-profile}
  \item{Imax}{The parameter given in the Morris down-profile}
  \item{design}{the matrix of the final design (low-discrepancy LHS)}
  \item{critValues}{vector of criterion values along the iterations}
  \item{tempValues}{vector of temperature values along the iterations}
  \item{probaValues}{vector of acceptation probability values along the iterations}
}

\references{

Damblin G., Couplet M., and Iooss B. (2013). Numerical studies of space filling designs: optimization of Latin Hypercube Samples and subprojection properties, \emph{Journal of Simulation,} 7:276-289, 2013.

M. Morris and J. Mitchell (1995) Exploratory designs for computationnal experiments. Journal of 
Statistical Planning and Inference, 43:381-402.

R. Jin, W. Chen and A. Sudjianto (2005) An efficient algorithm for constructing optimal design
of computer experiments. Journal of Statistical Planning and Inference, 134:268-287.

}

\author{G. Damblin & B. Iooss}

\seealso{Latin Hypercube Sample(\code{\link{lhsDesign}}),discrepancy criteria(\code{\link{discrepancyCriteria}}), geometric criterion (\code{\link{mindist}}\code{\link{phiP}}), optimization (\code{\link{maximinSA_LHS}},\code{\link{maximinESE_LHS}} ,\code{\link{discrepESE_LHS}})}


\examples{
dimension <- 2
n <- 10
X <- lhsDesign(n, dimension)$design

## Optimize the LHS with C2 criterion
Xopt <- discrepSA_LHS(X, T0=10, c=0.99, it=2000, criterion="C2")
plot(Xopt$design)
plot(Xopt$critValues, type="l")

## Optimize the LHS with C2 criterion and GEOM_MORRIS profile
\dontrun{
Xopt2 <- discrepSA_LHS(X, T0=10, c=0.99, it=1000, criterion="C2", profile="GEOM_MORRIS")
plot(Xopt2$design)
}
}

\keyword{design}

