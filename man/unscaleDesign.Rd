\name{unscaleDesign}
\alias{unscaleDesign}

\title{Unscale a Design}

\description{
     This function unscales the values of a scaled design (values in [0,1]).
     The unscaling can be made by the inverse Rosenblatt transformation (by applying the empirical quantile function given by another design) or by translating the design from maximum and minimum values (given for each variable).
}

\usage{
unscaleDesign(design, min=NULL, max=NULL, uniformize=FALSE, InitialDesign=NULL)
}

\arguments{
  \item{design}{a matrix (or a data.frame) corresponding to the design of experiments to unscale}
  \item{min}{the vector of minimal bounds of each design variable}
  \item{max}{the vector of maximal bounds of each design variable}
  \item{uniformize}{boolean: TRUE to use the inverse Rosenblatt transformation (the min and max vectors are useless in this case). If FALSE (default value), the translation from max and min values is applied}
  \item{InitialDesign}{If the inverse Rosenblatt transformation is applied (uniformize = TRUE): a matrix (or a data.frame) corresponding to the design which gives the empirical quantiles}
}

\value{A list containing:
  \item{design}{the unscaled design}
  \item{min}{the vector of minimal bounds that has been used}
  \item{max}{the vector of maximal bounds that has been used}
  \item{uniformize}{the value of this boolean argument}
  }

\author{B. Iooss}

\examples{
  d <- 2
  n <- 100
  x <- matrix(rnorm(d*n), ncol=d)
  xscale <- scaleDesign(x, uniformize=TRUE)
  xunscale1 <- unscaleDesign(xscale$design, uniformize=TRUE, InitialDesign=x)
  xunscale2 <- unscaleDesign(xscale$design, 
  min=c(min(x[,1]), min(x[,2])), max = c(max(x[,1]), max(x[,2])))
  par(mfrow=c(2,2))
  plot(x) ; plot(xscale$design)
  plot(xunscale1$design) ; plot(xunscale2$design) 
}
