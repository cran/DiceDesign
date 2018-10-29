olhDesign <- function(dimension, range = c(0, 1)) {
    funM <- function(dimpow2) {
        if (dimpow2 == 1) { matrix(c(1,2,2,1), 2, 2)
        } else {
            M <- funM(dimpow2 - 1)
            rbind(cbind(M, M+ncol(M)), cbind(M+ncol(M), M))
        }
    }
    funS <- function(dimpow2) {
        if (dimpow2 == 1) { matrix(c(1,1,1,-1), 2, 2)
        } else {
            S <- funS(dimpow2 - 1)
            P <- S[seq(1, nrow(S)/2),, drop = FALSE]
            Q <- S[seq(1+nrow(S)/2, nrow(S)),, drop = FALSE]
            rbind(cbind(P,P), cbind(Q,-Q), cbind(P,-P), cbind(Q,Q))
        }
    }
    funOLH <- function(dimpow2) {
        T <- funM(dimpow2) * funS(dimpow2)
        rbind(T, rep(0, ncol(T)), -T)
    }
    
    dimpow2 <- max(1, ceiling(log2(dimension)))
    mat <- funOLH(dimpow2)[, seq_len(dimension), drop = FALSE]
    mat <- if (identical(range, c(1,1))) mat
              else if (identical(range, c(0, 0))) mat +(nrow(mat)-1)/2  
                   else (mat/(nrow(mat)-1) +0.5) * diff(range) + range[1]
    dimnames(mat) <- NULL
    # dimnames(mat) <- list(seq_len(nrow(mat)), 
                          # paste0("X", seq_len(ncol(mat))))
list(n = nrow(mat), dimension = dimension, design = mat)
}

