
faureprimeDesign <- function(dimension, u = 2, range = c(0, -1)){
	if (!is.null(range)) {
		stopifnot(!any(c(NA, NaN, -Inf, Inf) %in% range))
		stopifnot(length(range) == 2L)
	}
	stopifnot(u >= 2)
	stopifnot(dimension >= 3)
	stopifnot(dimension <= 199)
	u         <- floor(u)
	dimension <- floor(dimension)
	primesupto199 <- c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,
					   71,73,79,83,89,97,101,103,107,109,113,127,131,137,
					   139,149,151,157,163,167,173,179,181,191,193,197,199)
	prime  <- primesupto199[primesupto199 >= dimension][1]
	stopifnot(prime^(u+1) <= 200^3)
	pu     <- prime^u
	n      <- pu-1
	design <- pu*runif.faure(n, dimension)$design
	if (!is.null(range)) {
		range <- as.numeric(range)
		if (identical(range, c(0,0))   | identical(range, c(0,1)) | 
			identical(range, c(0,-1))  | identical(range, c(1,1)) | 
			identical(range, c(-1,-1)) | identical(range, c(-1,1))) {
			# c(0, 0)  (0, n) already defined
			# c(1,1)   (1-n, n-1) = (2-p^u, p^u -2)
			# c(0,1)   (0, 1)
			# c(0,-1)  (p^-u, 1-p^-u)
			# c(-1,-1) (-1+2p^(-u) , 1-2p^(-u))
			# c(-1, 1) (-1, 1)
			if (identical(range, c(1,1)))   design <- 2*design -n-1  
			if (identical(range, c(0,1)))   design <- (-1 +design)/(n-1)
			if (identical(range, c(0,-1)))  design <- design/(n+1) 
			if (identical(range, c(-1,-1))) design <- (2*design -n-1)/(n+1)
			if (identical(range, c(-1,1)))  design <- (2*design -n-1)/(n-1)
		} else {
			stopifnot(range[2] > range[1])
			design <- (-1 +design)/(n-1)*diff(range) +range[1]
		}
	}
	list(design = design, n = (prime^u)-1, dimension = dimension, 
         prime = prime, u = u)
}

