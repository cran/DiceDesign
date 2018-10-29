nolhDesign <- function(dimension, range = c(0, 1)) {
    # NOLHdes <- get(load(system.file("data", "NOLHdesigns.rda", package = "DiceDesign")))
    NOLHdes   <- get("NOLHdesigns", pos = "package:DiceDesign")
    dfr <-  if (dimension %in% 2:7) NOLHdes[["nolh2_7"]]
                else if (dimension %in% 8:11) NOLHdes[["nolh8_11"]]
                    else if (dimension %in% 12:16) NOLHdes[["nolh12_16"]]
                        else if (dimension %in% 17:22) NOLHdes[["nolh17_22"]]
                            else if (dimension %in% 23:29) NOLHdes[["nolh23_29"]]
                                else stop("dimension must be in 2:29")
    # dfr <-  if (dimension %in% 2:7) NOLHdesigns[["nolh2_7"]]
                # else if (dimension %in% 8:11) NOLHdesigns[["nolh8_11"]]
                    # else if (dimension %in% 12:16) NOLHdesigns[["nolh12_16"]]
                        # else if (dimension %in% 17:22) NOLHdesigns[["nolh17_22"]]
                            # else if (dimension %in% 23:29) NOLHdesigns[["nolh23_29"]]
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

