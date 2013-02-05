\name{discrepESE_LHS}
\alias{discrepESE_LHS}

\title{Enhanced Stochastic Evolutionnary (ESE) algorithm for Latin Hypercube Sample (LHS) optimization via L2-discrepancy criteria}
\description{The objective is to produce low-discrepancy LHS. ESE is a powerful genetic algorithm to produce space-filling designs. It has been adapted here to main discrepancy criteria.}
\usage{discrepESE_LHS(design,T0=0.005*discrepancyCriteria(design)[[1]],
inner_it=100,J=50,it=2,criterion="DC2")}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments.}
  \item{T0}{The initial temperature.}
  \item{inner_it}{The number of iterations for inner loop}
  \item{J}{The number of new proposed LHS inside the inner loop}
  \item{it}{ The number of iterations for outer loop}
  \item{criterion}{The criterion to be optimized ("DL2","DW2","DC2"). "DL2" is the L2-star discrepancy ; "DW2" is the wrap-around L2-discrepancy ; DC2" is centered L2-discrepancy and is the choice by default}
  }

\details{This function implements a stochastic algorithm (ESE) to produce optimized LHS. It is based on Jin et al works (2005). Here, It has been adapted to some discrepancy criteria taking into account new ideas about the revaluations of discrepancy value after a LHS elementary perturbation (in order to avoid computing all terms in the discrepancy formulas).}


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

\seealso{discrepancy criteria(\code{\link{discrepancyCriteria}}), geometric criterion (\code{\link{mindist}}\code{\link{phiP}}), optimization (\code{\link{maximinSA_LHS}},\code{\link{maximinESE_LHS}} ,\code{\link{discrepSA_LHS}})}


\examples{
\dontrun{
dimension <- 2
n <- 10
library(lhs)
X <- randomLHS(n,dimension)
discrepESE_LHS(X,T0=0.005*discrepancyCriteria(X)[[1]],inner_it=100,J=50,it=2)
}
}

\keyword{ design }