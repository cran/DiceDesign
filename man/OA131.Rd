\name{OA131}
\alias{OA131}
\docType{data}
\title{ A 3D orthogonal array of strength 2}
\description{
  A 3-dimensional linear orthogonal array (OA) of strength 2 with 49 points. The design points are equally spaced into 2 dimensional coordinate planes. However by construction, such OAs satisfy a linear relation, here: \code{x1 + 3*x2 + x3 = 0 (mod 7)}. As a consequence, the design points are contained in parallel planes orthogonal to (1,3,1). Actually, they are also contained in parallel planes orthogonal to other directions, as (2,-1,2) or (3,2,3), since the congruence relation leads to \code{2*x1 - x2 + 2*x3 = 0 (mod 7)} or \code{3*x1 + 2*x2 + 3*x3 = 0 (mod 7)}. For instance, they are contained in 4 parallel planes orthogonal to (2,-1,2).
}

\usage{data(OA131)}
\format{
  A data frame with 49 observations on the following 3 variables.
  \describe{
    \item{\code{x1}}{first coordinate}
    \item{\code{x2}}{second coordinate}
    \item{\code{x3}}{third coordinate}
  }
}

\source{ 
Roustant O., Franco J., Carraro L., Jourdan A. (2010), A radial scanning statistic for selecting space-filling designs in computer experiments, MODA-9 proceedings, \url{http://www.emse.fr/~roustant/index.html}.
}

\examples{
data(OA131)

# centering and reducing to [0,1]^3
OA <- (OA131 + 0.5)/7
pairs(OA, xlim=c(0,1), ylim=c(0,1))

\dontrun{
library(lattice)
cloud(x3~x1+x2, data=OA, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), 
      screen = list(z = 50, x = -70, y = 0))}
}

\keyword{datasets}

