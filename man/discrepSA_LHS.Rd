\name{discrepSA_LHS}
\alias{discrepSA_LHS}

\title{Simulated annealing (SA) routine for Latin Hypercube Sample (LHS) optimization via L2-discrepancy criteria}
\description{The objective is to produce low-discrepancy LHS. SA is an efficient algorithm to produce space-filling designs. It has been adapted here to main discrepancy criteria.}
\usage{discrepSA_LHS(design,T0=10,c=0.95,it=2000,criterion="DC2",
profile="GEOM",Imax=100)}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments}
  \item{T0}{The initial temperature}
  \item{c}{A constant parameter regulating how the temperature goes down}
  \item{it}{The number of iterations}
  \item{criterion}{The criterion to be optimized. One can choose three different L2-discrepancies : the C2 (centered) discrepancy ("DC2"), the L2 star discrepancy ("DL2") and the W2 (wrap-around) discrepancy ("DW2")}
  \item{profile}{The temperature down-profile, geometric called "GEOM" or linear according to the Morris algorithm called "LINEAR_MORRIS"}
  \item{Imax}{A parameter given only if you choose the linear down-profile. It adjusts the number of iterations without improvement before a new elementary perturbation}
  }


\details{This function implements a classical routine to produce optimized LHS. It is based on the work of Morris and Mitchell (1995). They have proposed a SA version for LHS optimization according to mindist criterion. Here, it has been adapted to some discrepancy criteria taking in account new ideas about the reevaluations of a discrepancy value after a LHS elementary perturbation (in order to avoid computing all terms in the discrepancy formulas).}


\value{A low-discrepancy LHS design.}

\references{

Damblin G., Couplet M., and Iooss B. (2013). Numerical studies of space filling designs : optimization algorithms and subprojection properties, \emph{Journal of Simulation, submitted}. 
\url{http://www.gdr-mascotnum.fr/doku.php?id=iooss1}

D.Morris and J.Mitchell (1995) Exploratory designs for computationnal experiments. Journal of 
Statistical Planning and Inference, 43:381-402.

R.Jin, W.Chen, A.Sudjianto (2005) An efficient algorithm for constructing optimal design
of computer experiments. Journal of Statistical Planning and Inference, 134:268-287.

A.Marrel (2008) Mise en oeuvre et utilisation du meta-modele processus gaussien pour
l'analyse de sensibilite de modeles numeriques. PhD thesis, INSA Toulouse, France.

}

\author{G.Damblin & B. Iooss}

\seealso{discrepancy criteria(\code{\link{discrepancyCriteria}}), geometric criterion (\code{\link{mindist}}\code{\link{phiP}}), optimization (\code{\link{maximinSA_LHS}},\code{\link{maximinESE_LHS}} ,\code{\link{discrepESE_LHS}})}


\examples{
dimension <- 2
n <- 10
library(lhs)
X <- randomLHS(n,dimension)
discrepSA_LHS(X,T0=10,c=0.95,it=2000,criterion="DC2")
\dontrun{
  discrepSA_LHS(X,T0=10,c=0.95,it=2000,criterion="DC2",profile="LINEAR_MORRIS")
}
}

\keyword{ design }