\name{OA131_scrambled}
\alias{OA131_scrambled}
\docType{data}
\title{ A scrambled 3D orthogonal array of strength 2}
\description{
This design is obtained by adding a uniform noise to each coordinate of the orthogonal array OA131.
}

\usage{data(OA131_scrambled)}
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
data(OA131_scrambled)

pairs(OA131, xlim=c(0,1), ylim=c(0,1))
pairs(OA131_scrambled, xlim=c(0,1), ylim=c(0,1))
}

\keyword{datasets}

