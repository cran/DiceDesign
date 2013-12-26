\name{mstCriteria}
\alias{mstCriteria}

\title{Deriving the MST criteria}
\description{Compute both the mean and the standard deviation of the Minimal Spanning Tree (MST)}
\usage{mstCriteria(design,plot2d="FALSE")}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments.}
  \item{plot2d}{an argument for visualizing the mst of a 2d design} 
  }


\details{In our context, a MST is a tree whose the sum of the lengthes of the edges is minimal. Even if unicity does not hold, the overall length is stable. The mean and the standard deviation of the lengthes of the edges are usually derived to analyze the geometric profile of the design. A large mean and a small standard deviation characterize a so-called quasi-periodic design.}


\value{A list containing two components:
\item{tree}{a list containing the MST: each component of it contains a vector with all vertices which are connected with the experiment corresponding to the number of the components}
\item{stats}{vector with both the mean and the standard deviation values of the lengthes of the edges}
}

\references{

Damblin G., Couplet M., and Iooss B. (2013). Numerical studies of space filling designs: optimization of Latin hypercube samples and subprojection properties,
\emph{Journal of Simulation,} 7:276-289, 2013.
\url{http://www.gdr-mascotnum.fr/doku.php?id=iooss1}

Dussert, C., Rasigni, G., Rasigni, M., and Palmari, J. (1986). Minimal spanning tree: A new approach for studying order and disorder. Physical Review B, 34(5):3528-3531.

Franco J. (2008). Planification d'experiences numerique en phase exploratoire pour la simulation des phenomenes complexes, \emph{PhD thesis, Ecole Nationale Superieure des Mines de Saint Etienne}.

Franco, J., Vasseur, O., Corre, B., and Sergent, M. (2009). Minimum spanning tree: A new approach to assess the quality of the design of computer experiments. Chemometrics and Intelligent Laboratory Systems, 97:164-169.

Prim, R.C. (1957). Shortest connection networks and some generalizations, \emph{in Bell System Technical Journal 36:1389-1401}.

}

\author{G.Damblin & B.Iooss}


\examples{
dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
mstCriteria(X,plot2d=TRUE)
}

\keyword{ design }