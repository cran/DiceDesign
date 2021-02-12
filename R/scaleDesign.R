
scaleDesign <- function(design, min=NULL, max=NULL, uniformize=FALSE){
  d <- dim(design)[[2]] ; n <- dim(design)[[1]]
  
  if (uniformize){ # Rosenblatt transformation to U[0,1]^d
    x <- design
    for (i in 1:d) x[,i] <- ecdf(design[,i])(design[,i])
    
  }else{ # scaling between min and max
    if (is.null(min)) min <- apply(design, 2, min)
    if (is.null(max)) max <- apply(design, 2, max)
    x <- ( design - matrix(min, ncol=d, nrow=n, byrow=T) ) / matrix((max - min), ncol=d,nrow=n, byrow=T)
  }
  return(list(design = x, min = min, max = max, uniformize = uniformize, InitialDesign=design))
}
