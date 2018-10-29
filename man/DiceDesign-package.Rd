\name{DiceDesign-package}
\alias{DiceDesign-package}
\alias{DiceDesign}
\docType{package}
\title{Designs of Computer Experiments}
\description{Space-Filling Designs (SFD) and space-filling criteria (distance-based and uniformity-based).}
\details{
%\tabular{ll}{
%	Package: \tab DiceDesign\cr
%	Type: \tab Package\cr
%	Version: \tab 1.7\cr
%	Date: \tab 2015-06-15\cr
%	License: \tab GPL-3\cr
%}
This package provides tools to create some specific Space-Filling Design (SFD) and to test their quality:
\itemize{
    \item Latin Hypercube designs (randomized or centered)
    \item Strauss SFD and Maximum entropy SFD, WSP designs
    \item Optimal (low-discrepancy and maximin) Latin Hypercube desigsn by simulated annealing and genetic algorithms,
    \item Orthogonal and Nearly Orthogonal Latin Hypercube designs,
    \item Discrepancies criteria, distance measures,
    \item Minimal spanning tree criteria,
    \item Radial scanning statistic
}
}

\author{
J. Franco, D. Dupuy, O. Roustant, P. Kiener, G. Damblin and B. Iooss. Thanks to A. Jourdan for discussions about OA131.

Maintainer: Celine Helbert \email{Celine.Helbert@ec-lyon.fr} 
}

\note{
Part of this work was conducted within the frame of the DICE (Deep Inside Computer Experiments) Consortium between ARMINES, Renault, EDF, IRSN, ONERA and TOTAL S.A. (http://dice.emse.fr/).

In this package, only Faure's sequence is implemented. Note that the \pkg{randtoolbox} package provides the following quasi random sequences: the Sobol sequence, the Halton (hence Van Der Corput) sequence and the Torus sequence (also known as Kronecker sequence). Note also that the \pkg{lhs} package provides other types of algorithms to compute optimized LHS.
}

\references{
Cioppa T.M., Lucas T.W. (2007). Efficient nearly orthogonal and space-filling Latin hypercubes. \emph{Technometrics} 49, 45-55. \url{http://www.dtic.mil/dtic/tr/fulltext/u2/a520796.pdf}.

Damblin G., Couplet M., and Iooss B. (2013). Numerical studies of space filling designs: optimization of Latin Hypercube Samples and subprojection properties, \emph{Journal of Simulation}, 7:276-289, 2013. \url{http://www.gdr-mascotnum.fr/doku.php?id=iooss1}.

De Rainville F.-M., Gagne C., Teytaud O., Laurendeau D. (2012). Evolutionary optimization of low-discrepancy sequences. \emph{ACM Transactions on Modeling and Computer Simulation (TOMACS)}, 22(2), 9. \url{https://dl.acm.org/citation.cfm?id=2133393}.

Dupuy D., Helbert C., Franco J. (2015), DiceDesign and DiceEval: Two R-Packages for Design and Analysis of Computer Experiments, \emph{Journal of Statistical Software}, \bold{65}(11), 1--38, \url{http://www.jstatsoft.org/v65/i11/}.

Fang K.-T., Li R. and Sudjianto A. (2006) Design and Modeling for Computer Experiments, \emph{Chapman & Hall}.

Nguyen N.K. (2008) A new class of orthogonal Latinhypercubes, \emph{Statistics and Applications}, Volume 6, issues 1 and 2, pp.119-123.

Roustant O., Franco J., Carraro L., Jourdan A. (2010), A radial scanning statistic for selecting space-filling designs in computer experiments, \emph{MODA-9 proceedings}, \url{http://www.emse.fr/~roustant/index.html}.

Santner T.J., Williams B.J. and Notz W.I. (2003) The Design and Analysis of Computer Experiments, \emph{Springer}, 121-161.
}

\examples{
# **********************
# Designs of experiments
# **********************

# A maximum entropy design with 20 points in [0,1]^2
p <- dmaxDesign(20,2,0.9,200)
plot(p$design,xlim=c(0,1),ylim=c(0,1))

# Change the dimnames, adjust to range (-10, 10) and round to 2 digits
xDRDN(p, letter = "T", dgts = 2, range = c(-10, 10))

# ************************
# Criteria: L2-discrepancy
# ************************
dp <- discrepancyCriteria(p$design,type=c('L2','C2'))
# Coverage measure
covp <- coverage(p$design)

# *******************************
# Criteria: Minimal Spanning Tree
# *******************************
mstCriteria(p$design,plot2d=TRUE)

# ****************************************************************
# Radial scanning statistic: Detection of defects of Sobol designs
# ****************************************************************

# requires randtoolbox package
library(randtoolbox)

# in 2D
rss <- rss2d(design=sobol(n=20, dim=2), lower=c(0,0), upper=c(1,1),
	type="l", col="red")

# in 8D. All pairs of dimensions are tried to detect the worst defect
# (according to the specified goodness-of-fit statistic).
d <- 8
n <- 10*d
rss <- rss2d(design=sobol(n=n, dim=d), lower=rep(0,d), upper=rep(1,d),
	type="l", col="red")

# avoid this defect with scrambling ?
#    1. Faure-Tezuka scrambling (type "?sobol" for more details and options)
rss <- rss2d(design=sobol(n=n, dim=d, scrambling=2), lower=rep(0,d),
	upper=rep(1,d), type="l", col="red")
#    2. Owen scrambling
rss <- rss2d(design=sobol(n=n, dim=d, scrambling=1), lower=rep(0,d),
	upper=rep(1,d), type="l", col="red")

}
\keyword{package}

