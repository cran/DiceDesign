\name{DiceDesign-package}
\alias{DiceDesign-package}
\alias{DiceDesign}
\docType{package}
\title{Designs of Computer Experiments}
\description{Space-Filling Designs (SFD) and criteria of uniformity.}
\details{\tabular{ll}{
	Package: \tab DiceDesign\cr
	Type: \tab Package\cr
	Version: \tab 1.0\cr
	Date: \tab 2010-02-17\cr
	License: \tab GPL-3\cr
	}
This package provides tools to create some specific Space-Filling Design (SFD) and to test their quality:
\itemize{
    \item Strauss SFD and Maximum entropy SFD,
    \item Discrepancies criteria, distance measures,
    \item Radial scanning statistic}
}

\author{J. Franco, D. Dupuy and O. Roustant. Thanks to A. Jourdan for discussions about OA131.

(maintainer: D. Dupuy \email{dupuy@emse.fr})}

\note{This work was conducted within the frame of the DICE (Deep Inside Computer Experiments) Consortium between ARMINES, Renault, EDF, IRSN, ONERA and TOTAL S.A. (http://www.dice-consortium.fr/).

In this package only Faure's sequence is implemented.
Note that the \pkg{randtoolbox} package provides the following quasi
random sequences: the Sobol sequence, the Halton (hence Van
Der Corput) sequence and the Torus sequence (also known as
Kronecker sequence).

We refer to the \pkg{lhs} package to construct Latin Hypercube Designs.}

\references{
Fang K.-T., Li R. and Sudjianto A. (2006) Design and Modeling for Computer Experiments, \emph{Chapman & Hall}.

Santner T.J., Williams B.J. and Notz W.I. (2003) The Design and Analysis of Computer Experiments, \emph{Springer}, 121-161.

Roustant O., Franco J., Carraro L., Jourdan A. (2010), A radial scanning statistic for selecting space-filling designs in computer experiments, MODA-9 proceedings, \url{http://www.emse.fr/~roustant}

}

\examples{
# **********************
# Designs of experiments
# **********************

# A maximum entropy design with 20 points in [0,1]^2
p <- dmaxDesign(20,2,0.9,200)
plot(p$design,xlim=c(0,1),ylim=c(0,1))

# ************************
# Criteria: L2-discrepancy
# ************************
dp <- discrepancyCriteria(p$design,type=c('L2','C2'))
# Coverage measure
covp <- coverage(p$design)

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
\keyword{ package }