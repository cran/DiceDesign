nolhdrDesign <- function(dimension, range = c(0, 1)) {
    # NOLHDRdes <- get(load(system.file("data", "NOLHDRdesigns.rda", package = "DiceDesign")))
    NOLHdes   <- get("NOLHdesigns", pos = "package:DiceDesign")
    NOLHDRdes <- get("NOLHDRdesigns", pos = "package:DiceDesign")
    dfr <-  if (dimension %in% 2:7) NOLHdes[["nolh2_7"]]
                else if (dimension %in% 8:29) NOLHDRdes[[dimension - 7]]
                    else stop("dimension must be in 2:29")
    # dfr <-  if (dimension %in% 2:7) NOLHdesigns[["nolh2_7"]]
                # else if (dimension %in% 8:29) NOLHDRdesigns[[dimension - 7]]
                    # else stop("dimension must be in 2:29")
    mat <- as.matrix(dfr[,seq_len(dimension)])
    mat <- if (identical(range, c(1,1))) mat 
              else if (identical(range, c(0, 0))) mat +(nrow(mat)-1)/2  
                   else (mat/(nrow(mat)-1) +0.5) *diff(range) + range[1]
    dimnames(mat) <- NULL
    # dimnames(mat) <- list(seq_len(nrow(mat)), 
                          # paste0("X", seq_len(ncol(mat))))
list(n = nrow(mat), dimension = dimension, design = mat)
}

