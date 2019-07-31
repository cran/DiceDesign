\name{scaleDesign}
\alias{scaleDesign}

\title{Scale a Design}

\description{
     This function scales the values of the design points to values comprised in [0,1].
     The scaling can be made by the Rosenblatt transformation (uniformization by applying the empirical cumulative distribution function) or by translating the design from maximum and minimum values (given for each variable).
}

\usage{
scaleDesign(design, min=NULL, max=NULL, uniformize=FALSE)
}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments to scale}
  \item{min}{the vector of minimal bounds of each design variable. If not given, the minimal value of each variable is taken}
  \item{max}{the vector of maximal bounds of each design variable. If not given, the maximal value of each variable is taken}
  \item{uniformize}{boolean: TRUE to use the Rosenblatt transformation (the min and max vectors are useless in this case). If FALSE (default value), the translation from max and min values is applied}
}

\value{A list containing:
  \item{design}{the scaled design}
  \item{min}{the vector of minimal bounds that has been used}
  \item{max}{the vector of maximal bounds that has been used}
  \item{uniformize}{the value of this boolean argument}
  \item{InitialDesign}{the starting design}
  }

\author{B. Iooss}

\examples{
  d <- 2
  n <- 100
  x <- matrix(rnorm(d*n), ncol=d)
  xscale1 <- scaleDesign(x, uniformize=FALSE)
  xscale2 <- scaleDesign(x, uniformize=TRUE)
  par(mfrow=c(1,2))
  plot(xscale1$design) ; plot(xscale2$design)
}
