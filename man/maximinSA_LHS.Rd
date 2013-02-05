\name{maximinSA_LHS}
\alias{maximinSA_LHS}

\title{Simulated annealing (SA) routine for Latin Hypercube Sample (LHS) optimization via phiP criteria}
\description{The objective is to produce maximin LHS. SA is an efficient algorithm to produce space-filling designs.}
\usage{maximinSA_LHS(design,T0=10,c=0.95,it=2000,p=50,profile="GEOM",Imax=100)}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments}
  \item{T0}{The initial temperature}
  \item{c}{A constant parameter regulating how the temperature goes down}
  \item{it}{The number of iterations}
  \item{p}{power required in phiP criterion}
  \item{profile}{The temperature down-profile, geometric called "GEOM" or linear according to the Morris algorithm called "LINEAR_MORRIS"}
  \item{Imax}{A parameter given only if you choose the linear down-profile. It adjusts the number of iterations without improvement before a new elementary perturbation}
  }


\details{This function implements a classical routine to produce optimized LHS. It is based on the work of Morris and Mitchell (1995). They have proposed a SA version for LHS optimization according to mindist criterion. Here, it has been adapted to the \code{phiP} criterion. It has been shown (Pronzato and Muller, 2012, Damblin et al., 2013) that optimizing phiP is more efficient to produce maximin designs than optimizing mindist. When \eqn{p} tends to infinity, optimizing a design with \code{phi_p} is equivalent to optimizing a design with \code{mindist}.}


\value{A maximin LHS design.}

\references{

Damblin G., Couplet M., and Iooss B. (2013). Numerical studies of space filling designs : optimization algorithms and subprojection properties, \emph{Journal of Simulation, submitted}. 
\url{http://www.gdr-mascotnum.fr/doku.php?id=iooss1}

D.Morris and J.Mitchell (1995) Exploratory designs for computationnal experiments. Journal of 
Statistical Planning and Inference, 43:381-402.

R.Jin, W.Chen, A.Sudjianto (2005) An efficient algorithm for constructing optimal design
of computer experiments. Journal of Statistical Planning and Inference, 134:268-287.

A.Marrel (2008) Mise en oeuvre et utilisation du meta-modele processus gaussien pour
l'analyse de sensibilite de modeles numeriques. PhD thesis, INSA Toulouse, France.

Pronzato, L. and Muller, W. (2012). Design of computer experiments: space filling and beyond, \emph{Statistics and Computing,} 22:681-701.
}

\author{G.Damblin & B. Iooss}

\seealso{discrepancy criteria(\code{\link{discrepancyCriteria}}), geometric criterion (\code{\link{mindist}}\code{\link{phiP}}), optimization (\code{\link{discrepSA_LHS}},\code{\link{maximinESE_LHS}} ,\code{\link{discrepESE_LHS}})}


\examples{
dimension <- 2
n <- 10
library(lhs)
X <- randomLHS(n,dimension)
maximinSA_LHS(X,T0=10,c=0.95,it=2000)
\dontrun{
  maximinSA_LHS(X,T0=10,c=0.95,it=2000,profile="LINEAR_MORRIS")
}
}

\keyword{ design }