\name{straussDesign}
\alias{straussDesign}

\title{ Designs based on Strauss process}
\description{Space-Filling Designs based on Strauss process}
\usage{
straussDesign(n,dimension, RND, alpha=0.5, repulsion=0.001, NMC=1000,
              constraints1D=0, repulsion1D=0.0001, seed=NULL)
}

\arguments{
  \item{n        }{the number of experiments}
  \item{dimension}{the number of input variables}
  \item{RND      }{a real number which represents the radius of interaction} 
  \item{alpha    }{the potential power (default, fixed at 0.5)}
  \item{repulsion}{the repulsion parameter in the unit cube (gamma)}
  \item{NMC      }{the number of McMC iterations (this number must be large to converge)}
  \item{constraints1D	}{1 to impose 1D projection constraints, 0 otherwise}
  \item{repulsion1D }{the repulsion parameter in 1D}
  \item{seed	 }{ seed for the uniform generation of number}}

\details{Strauss designs are Space-Filling designs initially defined from Strauss process:
\deqn{\pi (X) = k \gamma^{s(X)}}{\pi (X) = k gamma^s(X)}
where \eqn{s(X)} is is the number of pairs of points \eqn{(x^{i}, x^{j})}{(xi,xj)} of the design \eqn{X = \left( x^{1}, \ldots, x^{n} \right)}{X = ( x1, ..., xn \right)} that are separated by a distance no greater than the radius of interaction \code{RND}, \eqn{k} is the normalizing constant and \eqn{\gamma}{gamma} is the repulsion parameter. This distribution corresponds to the particular case \code{alpha}=0.

For the general case, a stochastic simulation is used to construct a Markov chain which converges to a spatial density of points \eqn{\pi(X)}{\pi(X)} described by the Strauss-Gibbs potential. In practice, the Metropolis-Hastings algorithm is implemented to simulate a distribution of points which converges to the stationary law: 
\deqn{\pi (X) \propto exp(-U(X))}{\pi(X) = k exp(-U(X))}
with a potentiel \eqn{U} defined by:
\deqn{U(X) = \beta \sum_{1 \leq i < j \leq n} \varphi  \left( \| x^{i}-x^{j} \| \right)}{beta Sum_{i<j} phi(|| xi-xj ||)}
where \eqn{\beta = - \ln \gamma}{beta = -ln(gamma)}, \eqn{\varphi (h) = \left( 1 - \frac{h}{RND} \right) ^{\alpha}}{phi (h) = (1-h/RND)^{alpha}} if \eqn{h \leq}{h <=} \code{RND} and 0 otherwise.


The input parameters of \code{straussDesign} function can be interpreted as follows: 

  - \code{RND} is used to compute the number of pairs of points of the design separated by a distance no more than \code{RND}. A point is said "in interaction" with another if the spheres of radius \code{RND}/2 centered on these points intersect.

  - \code{alpha} is the potential power \eqn{\alpha}{alpha}. The case \code{alpha}=0 corresponds to Strauss process (0-1 potential). 

  - \code{repulsion} is equal to the \eqn{\gamma}{gamma} parameter of the Strauss process. Note that \eqn{\gamma}{gamma} belongs to ]0,1].

  - \code{constraints1D} allows to specify some constraints into the margin. If \code{constraints1D}==1, two repulsion parameters are needed: one for the all space (\code{repulsion}) and the other for the 1D projection (\code{repulsion1D}). Default values are \code{repulsion}=0.001 and \code{repulsion1D}=0.001. Note that the value of the radius of interaction in the one-dimensional axis is not an input parameter and is automatically fixed at \eqn{0.75/n}.}

\value{A list containing:
  \item{n         }{the number of experiments}
  \item{dimension }{the number \eqn{d} of variables}
  \item{design_init}{the initial distribution of \eqn{n} points \eqn{[0,1]^{d}}}
  \item{radius    }{the radius of interaction}
  \item{alpha	  }{the potential power alpha}
  \item{repulsion }{the repulsion parameter \eqn{\gamma}{\gamma}}
  \item{NMC       }{the number of iterations McMC}
  \item{constraints1D  }{an integer indicating if constraints on the factorial axis are imposed. If its value is different from zero, a component \code{repulsion1D} containing the value of the repulsion parameter \eqn{\gamma}{gamma} in dimension 1 is added at the list.}
  \item{design    }{the design of experiments in [0,1]\eqn{^{d}}}
  \item{seed	  }{the seed corresponding to the design}
}
\author{J. Franco}

\references{
J. Franco, X. Bay, B. Corre and D. Dupuy (2008) Planification d'experiences numeriques a partir du processus ponctuel de Strauss, \url{https://hal.archives-ouvertes.fr/hal-00260701/fr/}.
}

\examples{
## Strauss-Gibbs designs in dimension 2 (n=20 points)
S1 <- straussDesign(n=20, dimension=2, RND=0.2)
plot(S1$design, xlim=c(0,1), ylim=c(0,1))
theta <- seq(0,2*pi, by=2*pi/(100 - 1))
for(i in 1:S1$n){
   lines(S1$design[i,1]+S1$radius/2*cos(theta),
	     S1$design[i,2]+S1$radius/2*sin(theta), col='red')
}

## 2D-Strauss design
S2 <- straussDesign(n=20, dimension=2, RND=0.2, NMC=200,
	         constraints1D=0, alpha=0, repulsion=0.01)
plot(S2$design,xlim=c(0,1),ylim=c(0,1))

## 2D-Strauss designs with constraints on the axis
S3 <- straussDesign(n=20, dimension=2, RND=0.18, NMC=200, 
	         constraints1D=1, alpha=0.5, repulsion=0.1, repulsion1D=0.01)
plot(S3$design, xlim=c(0,1),ylim=c(0,1))
rug(S3$design[,1], side=1)
rug(S3$design[,2], side=2)

## Change the dimnames, adjust to range (-10, 10) and round to 2 digits
xDRDN(S3, letter="T", dgts=2, range=c(-10, 10))
}

\keyword{design}

