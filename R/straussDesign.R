	straussDesign <- function(n,dimension,RND,NMC=1000,alpha=0.5,constraints1D=0,
			repulsion=10,gamma1D=0.1){
	 
	#-- Computation of maximal value for the radius in dimension ND
	#if (dimension%%2==0){
	#	k   <- floor(dimension/2);
	#	rnd <- (factorial(k) * 2^(dimension)/(n*pi^k))^(1/dimension);
	#}else{
	#	k   <- floor((dimension+1)/2);
	#	rnd <- (factorial(dimension) / (n*factorial(k-1)*pi^(k-1)))^(1/dimension);
	#}
	#if (RND>=rnd){
	#	stop(paste("The value of RND (here ",RND,") should be less than ",round(rnd,digits=4),".\n",sep=""),call.=TRUE); 
	#}
	
	cat("The radius of interactions in dimension ",dimension," is ",RND,"\n",sep="") 
	# Case with none constraints
	if(constraints1D==0){
		# Case with potential different to zero
		if (alpha != 0){
			betaND <- repulsion
			cat("Repulsion parameter (betaND) ",betaND,".\n",sep="")
		} else{ # Case with potential 0-1
			gammaND <- repulsion
			cat("Repulsion parameter (gammaND) ",gammaND,".\n",sep="")
		} # End of the case with none constraint
	}else{ # Constraint 1D as Hypercubes Latin
		R1D <- 0.75/n 
		# Case with potential different to zero
		if (alpha != 0){
			betaND <- repulsion
			cat("Repulsion parameter in dimension ",dimension," (Beta ND) = ",betaND,".\n",sep="")
			cat("Repulsion parameter for 1D constraints (Gamma 1D) = ",gamma1D,".\n",sep="")
			} else { # Case with potential 0-1
			gammaND <- repulsion
			cat("Repulsion parameter in dimension ",dimension," (Gamma ND) = ",gammaND,".\n",sep="")
			cat("Repulsion parameter for 1D constraints (Gamma 1D) = ",gamma1D,".\n",sep="")
		}
	}
		
	# intial design
	init <- runif(n*dimension)

	# plot(t(matrix(init,ncol=n, nrow=dimension,byrow=FALSE)),xlab='x1',ylab='x2',xlim=c(0,1),ylim=c(0,1),
	#		col='red',pch=20,xaxs = "i",yaxs = "i")
	
 	out <- .C("C_StraussDesign",
	 	as.double(init),as.integer(n),as.integer(dimension),as.integer(constraints1D),
	 	as.integer(NMC),as.double(RND),as.double(alpha),as.double(repulsion),
		as.double(gamma1D),ans = double(n * dimension),PACKAGE="DiceDesign")
	
  # Outputs
    return(list(n=n,dimension=dimension,design=matrix(out$ans,nrow=n, ncol=dimension,byrow=TRUE),
		radius=RND,NMC=NMC))
}