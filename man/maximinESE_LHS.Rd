\name{maximinESE_LHS}
\alias{maximinESE_LHS}

\title{Enhanced Stochastic Evolutionnary (ESE) algorithm for Latin Hypercube Sample (LHS) optimization via phiP criteria}
\description{The objective is to produce maximin LHS. ESE is a powerful genetic algorithm allowing to produce space-filling designs.}
\usage{maximinESE_LHS(design,T0=0.005*phiP(design,p=50),inner_it=100,J=50,it=1,p=50)}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments.}
  \item{T0}{The initial temperature.}
  \item{inner_it}{The number of iterations for inner loop}
  \item{J}{The number of new proposed LHS inside the inner loop}
  \item{it}{ The number of iterations for outer loop}
  \item{p}{power required in phiP criterion}
  }


\details{This function implements a stochastic algorithm (ESE) to produce optimized LHS. It is based on Jin et al works (2005).}


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

\seealso{discrepancy criteria(\code{\link{discrepancyCriteria}}), geometric criterion (\code{\link{mindist}}\code{\link{phiP}}), optimization (\code{\link{maximinSA_LHS}},\code{\link{discrepESE_LHS}} ,\code{\link{discrepSA_LHS}})}


\examples{
dimension <- 2
n <- 10
library(lhs)
X <- randomLHS(n,dimension)
maximinESE_LHS(X,T0=0.005*phiP(X),inner_it=100,J=50,it=2)

}

\keyword{ design }