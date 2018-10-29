\name{xDRDN}
\alias{xDRDN}

\title{Extract a Design and Give it a Range and Dimnames}
\description{
     Extract a design contained in a list (i.e. with a \eqn{design} item), adjust the range, give it dimnames and finally round the values to a certain number of digits. Colnames will look like (\code{A,B,C}), (\code{X1,X2,X3}), (\code{X01,X02,X03}), (\code{X001,X002,X003}). 
}
\usage{
xDRDN(obj, width = 1, letter = "X", dgts = NULL, range = NULL)
}

\arguments{
  \item{obj   }{ a list that contains a \code{design} item. Matrix or data.frame are also accepted }
  \item{width }{ the digit width in colnames (to write for instance  X1, X01, X001). If 0, colnames are filled with capital and small letters (without letters \code{I} and \code{i}) up to 50 columns }
  \item{letter}{ the generic letter used in colnames }
  \item{dgts  }{ the number of digits to which the design is rounded }
  \item{range }{ a vector \code{c(min, max)} to adjust the range of the design. The default \code{NULL} keeps the original range. Special ranges \code{c(0, 0)} and \code{c(1, 1)} are not accepted }
}

\value{
A rounded matrix or a data.frame with appropriate dimnames and an adjusted range. 
}

\examples{
xDRDN(lhsDesign(5, 12))
xDRDN(lhsDesign(5, 12), width = 2, letter = "V", dgts = 2, range = c(-10, 10))
head(xDRDN(olhDesign(50, range = c(1,1)), width = 0, letter = "Z"), 3)
head(xDRDN(olhDesign(51, range = c(1,1)), width = 0, letter = "Z"), 3)
}

