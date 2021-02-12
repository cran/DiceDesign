
unscaleDesign <- function(design, min=NULL, max=NULL, uniformize=FALSE, InitialDesign=NULL){
  d <- dim(design)[[2]] ; n <- dim(design)[[1]]
  
  if (uniformize & !is.null(InitialDesign)){ # Inverse Rosenblatt transformation from U[0,1]^d
    for (i in 1:d) x[,i] <- quantile(InitialDesign[,i],design[,i])
  }else x <- design * matrix((max - min), ncol=d,nrow=n, byrow=T) + matrix(min, ncol=d, nrow=n,byrow=T)
  
  return(list(design = x, min = min, max = max, uniformize = uniformize))
}
