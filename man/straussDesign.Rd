\name{straussDesign}
\alias{straussDesign}

\title{ Designs based on Strauss process}
\description{Space-Filling Designs based on Strauss process}
\usage{straussDesign(n,dimension,RND,NMC=1000,alpha=0.5,
		constraints1D=0,repulsion=10,gamma1D=0.1)
}

\arguments{
  \item{n        }{the number of experiments}
  \item{dimension}{the number of input variables}
  \item{RND      }{a real number which represents the radius of interaction} 
  \item{NMC      }{the number of McMC iterations (this number must be large to converge)}
  \item{alpha    }{the potential power (default, fixed at 0.5)}
  \item{constraints1D	}{1 to impose 1D projection constraints, 0 otherwise}
  
  \item{repulsion}{the repulsion parameter in unit cube}
  \item{gamma1D  }{the repulsion parameter in 1D (value between 0 and 1)}
}

\details{Strauss designs are Space-Filling designs initially defined from Strauss process. Through Markov chain Monte-Carlo methods, more general Gibbs processes can be used to perform different goals (for example constraints over the 1D projections)

  \code{RND} is used to compute the number of pairs of points of the design separated by a distance no more than RND.

  \code{alpha} is a parameter of interaction potential. The case \code{alpha}=0 corresponds to Gibbs process (0-1 potential). 

  \code{repulsion} is equal to \eqn{gamma}{gamma}=0.1 if alpha = 0 and to \eqn{beta=-ln(gamma)=10}{beta=-ln(gamma)=10} else.

  \code{constraints1D} permits to specify some constraints into the margin. If \code{constraints1D}==1, two repulsion parameters are needed: one for the all space (\code{repulsion}) and the other for the 1D projection (gamma1D). Default values are \code{betaND}=0.01 and \code{gamma1D}=0.1}

\value{A list containing:
  \item{n         }{the number of experiments}
  \item{dimension }{the number of variables}
  \item{design    }{the design of experiments in [0,1]\eqn{^{dimension}}{^{dimension}}}
  \item{radius    }{the radius of interaction}
  \item{NMC       }{the number of iterations McMC}
}
\author{J. Franco}

\references{J. Franco, X. Bay, B. Corre and D. Dupuy (2008) \emph{Planification
	d experiences numeriques a partir du processus ponctuel de Strauss}
	\url{http://hal.archives-ouvertes.fr/hal-00260701/fr/} (March 06,2008).}

\examples{
# Strauss-Gibbs designs in dimension 2 (n=20 points)
S1 <- straussDesign(n=20,dimension=2,RND=0.2)

plot(S1$design,xlim=c(0,1),ylim=c(0,1))
theta <- seq(0,2*pi,by =2*pi/(100 - 1))
for(i in 1:S1$n){
   lines(S1$design[i,1]+S1$radius/2*cos(theta),
	   S1$design[i,2]+S1$radius/2*sin(theta),col='red')
}
# 2D-Strauss design
S2 <- straussDesign(n=20,dimension=2,RND=0.2,NMC=200,
	constraints1D=0,alpha=0,repulsion=0.01)

plot(S2$design,xlim=c(0,1),ylim=c(0,1))

# 2D-Strauss designs with constraints on the axis
S3 <- straussDesign(n=20,dimension=2,RND=0.18,NMC=200,
	constraints1D=1,alpha=0.5,repulsion=10,gamma1D=0.1)

plot(S3$design,xlim=c(0,1),ylim=c(0,1))
rug(S3$design[,1],side=1)
rug(S3$design[,2],side=2)
}

\keyword{design }