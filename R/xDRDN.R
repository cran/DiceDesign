xDRDN <- function(obj, width = 1, letter = "X", dgts = NULL, range = NULL) {
    mat <-  if (is.list(obj) & "design" %in% names(obj)) {
                obj$design 
            } else { 
                if (is.matrix(obj) | is.data.frame(obj)) { 
                    obj 
                } else {
                    stop("obj must be a matrix, a data.frame or a list with a 'design' field") 
                }
            }
    style <- paste0("%", width, ".", width, "d")
    LET   <- c(LETTERS[c(1:8, 10:26)], letters[c(1:8, 10:26)])
    cn    <- if(width == 0 & ncol(mat) < 51) {
                   LET[seq_len(ncol(mat))]
             } else { 
                    paste0(letter, sprintf(style, seq_len(ncol(mat)))) 
             }
    dimnames(mat) <- list(seq_len(nrow(mat)), cn)
    currange <- range(mat)
    if (!is.null(range)) mat <- (mat - currange[1])/diff(currange)*diff(range) + range[1]
    if (!is.null(dgts))  mat <- round(mat, dgts)
return(mat)
}

