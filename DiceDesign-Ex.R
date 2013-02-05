pkgname <- "DiceDesign"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('DiceDesign')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("DiceDesign-package")
### * DiceDesign-package

flush(stderr()); flush(stdout())

### Name: DiceDesign-package
### Title: Designs of Computer Experiments
### Aliases: DiceDesign-package DiceDesign
### Keywords: package

### ** Examples

# **********************
# Designs of experiments
# **********************

# A maximum entropy design with 20 points in [0,1]^2
p <- dmaxDesign(20,2,0.9,200)
plot(p$design,xlim=c(0,1),ylim=c(0,1))

# ************************
# Criteria: L2-discrepancy
# ************************
dp <- discrepancyCriteria(p$design,type=c('L2','C2'))
# Coverage measure
covp <- coverage(p$design)

# *******************************
# Criteria: Minimal Spanning Tree
# *******************************
mstCriteria(p$design,plot2d=TRUE)

# ****************************************************************
# Radial scanning statistic: Detection of defects of Sobol designs
# ****************************************************************

# requires randtoolbox package
library(randtoolbox)

# in 2D
rss <- rss2d(design=sobol(n=20, dim=2), lower=c(0,0), upper=c(1,1),
	type="l", col="red")

# in 8D. All pairs of dimensions are tried to detect the worst defect
# (according to the specified goodness-of-fit statistic).
d <- 8
n <- 10*d
rss <- rss2d(design=sobol(n=n, dim=d), lower=rep(0,d), upper=rep(1,d),
	type="l", col="red")

# avoid this defect with scrambling ?
#    1. Faure-Tezuka scrambling (type "?sobol" for more details and options)
rss <- rss2d(design=sobol(n=n, dim=d, scrambling=2), lower=rep(0,d),
	upper=rep(1,d), type="l", col="red")
#    2. Owen scrambling
rss <- rss2d(design=sobol(n=n, dim=d, scrambling=1), lower=rep(0,d),
	upper=rep(1,d), type="l", col="red")




cleanEx()
nameEx("OA131")
### * OA131

flush(stderr()); flush(stdout())

### Name: OA131
### Title: A 3D orthogonal array of strength 2
### Aliases: OA131
### Keywords: datasets

### ** Examples

data(OA131)

# centering and reducing to [0,1]^3
OA <- (OA131 + 0.5)/7

pairs(OA, xlim=c(0,1), ylim=c(0,1))
## Not run: 
##D library(lattice)
##D cloud(x3~x1+x2, data=OA, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), 
##D     screen = list(z = 50, x = -70, y = 0))
## End(Not run)



cleanEx()
nameEx("OA131_scrambled")
### * OA131_scrambled

flush(stderr()); flush(stdout())

### Name: OA131_scrambled
### Title: A scrambled 3D orthogonal array of strength 2
### Aliases: OA131_scrambled
### Keywords: datasets

### ** Examples

data(OA131)
data(OA131_scrambled)

pairs(OA131, xlim=c(0,1), ylim=c(0,1))
pairs(OA131_scrambled, xlim=c(0,1), ylim=c(0,1))



cleanEx()
nameEx("coverage")
### * coverage

flush(stderr()); flush(stdout())

### Name: coverage
### Title: Coverage
### Aliases: coverage
### Keywords: design

### ** Examples

dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
coverage(X)



cleanEx()
nameEx("discrepESE_LHS")
### * discrepESE_LHS

flush(stderr()); flush(stdout())

### Name: discrepESE_LHS
### Title: Enhanced Stochastic Evolutionnary (ESE) algorithm for Latin
###   Hypercube Sample (LHS) optimization via L2-discrepancy criteria
### Aliases: discrepESE_LHS
### Keywords: design

### ** Examples

## Not run: 
##D dimension <- 2
##D n <- 10
##D library(lhs)
##D X <- randomLHS(n,dimension)
##D discrepESE_LHS(X,T0=0.005*discrepancyCriteria(X)[[1]],inner_it=100,J=50,it=2)
## End(Not run)



cleanEx()
nameEx("discrepSA_LHS")
### * discrepSA_LHS

flush(stderr()); flush(stdout())

### Name: discrepSA_LHS
### Title: Simulated annealing (SA) routine for Latin Hypercube Sample
###   (LHS) optimization via L2-discrepancy criteria
### Aliases: discrepSA_LHS
### Keywords: design

### ** Examples

dimension <- 2
n <- 10
library(lhs)
X <- randomLHS(n,dimension)
discrepSA_LHS(X,T0=10,c=0.95,it=2000,criterion="DC2")
## Not run: 
##D   discrepSA_LHS(X,T0=10,c=0.95,it=2000,criterion="DC2",profile="LINEAR_MORRIS")
## End(Not run)



cleanEx()
nameEx("discrepancyCriteria")
### * discrepancyCriteria

flush(stderr()); flush(stdout())

### Name: discrepancyCriteria
### Title: Discrepancy measure
### Aliases: discrepancyCriteria
### Keywords: design

### ** Examples

dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
discrepancyCriteria(X)



cleanEx()
nameEx("dmaxDesign")
### * dmaxDesign

flush(stderr()); flush(stdout())

### Name: dmaxDesign
### Title: Maximum Entropy Designs
### Aliases: dmaxDesign
### Keywords: design

### ** Examples

n <- 20
dimension <- 2
range <-0.9
niter_max <- 200
out <- dmaxDesign(n,dimension,range,niter_max)



cleanEx()
nameEx("factDesign")
### * factDesign

flush(stderr()); flush(stdout())

### Name: factDesign
### Title: Full Factorial Designs
### Aliases: factDesign
### Keywords: design

### ** Examples

## First example
g <- factDesign(2,7)
plot(g$design,xlim=c(0,1),ylim=c(0,1))
## Second example
g <- factDesign(2,c(2,7))
plot(g$design,xlim=c(0,1),ylim=c(0,1))



cleanEx()
nameEx("maximinESE_LHS")
### * maximinESE_LHS

flush(stderr()); flush(stdout())

### Name: maximinESE_LHS
### Title: Enhanced Stochastic Evolutionnary (ESE) algorithm for Latin
###   Hypercube Sample (LHS) optimization via phiP criteria
### Aliases: maximinESE_LHS
### Keywords: design

### ** Examples

dimension <- 2
n <- 10
library(lhs)
X <- randomLHS(n,dimension)
maximinESE_LHS(X,T0=0.005*phiP(X),inner_it=100,J=50,it=2)




cleanEx()
nameEx("maximinSA_LHS")
### * maximinSA_LHS

flush(stderr()); flush(stdout())

### Name: maximinSA_LHS
### Title: Simulated annealing (SA) routine for Latin Hypercube Sample
###   (LHS) optimization via phiP criteria
### Aliases: maximinSA_LHS
### Keywords: design

### ** Examples

dimension <- 2
n <- 10
library(lhs)
X <- randomLHS(n,dimension)
maximinSA_LHS(X,T0=10,c=0.95,it=2000)
## Not run: 
##D   maximinSA_LHS(X,T0=10,c=0.95,it=2000,profile="LINEAR_MORRIS")
## End(Not run)



cleanEx()
nameEx("meshRatio")
### * meshRatio

flush(stderr()); flush(stdout())

### Name: meshRatio
### Title: MeshRatio measure
### Aliases: meshRatio
### Keywords: design

### ** Examples

dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
meshRatio(X)



cleanEx()
nameEx("mindist")
### * mindist

flush(stderr()); flush(stdout())

### Name: mindist
### Title: Mindist measure
### Aliases: mindist
### Keywords: design

### ** Examples

dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
mindist(X)



cleanEx()
nameEx("mstCriteria")
### * mstCriteria

flush(stderr()); flush(stdout())

### Name: mstCriteria
### Title: Deriving the MST criteria
### Aliases: mstCriteria
### Keywords: design

### ** Examples

dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
mstCriteria(X,plot2d=TRUE)



cleanEx()
nameEx("phiP")
### * phiP

flush(stderr()); flush(stdout())

### Name: phiP
### Title: phiP criterion
### Aliases: phiP
### Keywords: design

### ** Examples

dimension <- 2
n <- 40
X <- matrix(runif(n*dimension),n,dimension)
phiP(X)



cleanEx()
nameEx("rss2d")
### * rss2d

flush(stderr()); flush(stdout())

### Name: rss2d
### Title: 2D graphical tool for defect detection of Space-Filling Designs.
### Aliases: rss2d
### Keywords: design

### ** Examples

# Detection of defects of Sobol designs

# requires randtoolbox package
library(randtoolbox)

	# in 2D
rss <- rss2d(design=sobol(n=20, dim=2), lower=c(0,0), upper=c(1,1), type="l", 
   col="red")

	# in 8D. 
	# All pairs of dimensions are tried to detect the worst defect 
	# (according to the specified goodness-of-fit statistic).
d <- 8
n <- 10*d
rss <- rss2d(design=sobol(n=n, dim=d), lower=rep(0,d), upper=rep(1,d), type="l", 
   col="red")

# avoid this defect with scrambling ?
#    1. Faure-Tezuka scrambling (type "?sobol" for more details and options)
rss <- rss2d(design=sobol(n=n, dim=d, scrambling=2), lower=rep(0,d), 
   upper=rep(1,d), type="l", col="red")
#    2. Owen scrambling
rss <- rss2d(design=sobol(n=n, dim=d, scrambling=1), lower=rep(0,d), 
   upper=rep(1,d), type="l", col="red")



cleanEx()
nameEx("rss3d")
### * rss3d

flush(stderr()); flush(stdout())

### Name: rss3d
### Title: 3D graphical tool for defect detection of Space-Filling Designs.
### Aliases: rss3d
### Keywords: design

### ** Examples

	# An orthogonal array in 3D   
data(OA131)

	# centering the design points of this 7-levels design
OA <- (OA131 + 0.5)/7

	# 2D projections onto coordinate axis 
pairs(OA, xlim=c(0,1), ylim=c(0,1))
	# Now let us look at the 3D properties with the 3D RSS (requires the rgl package)
rss <- rss3d(OA, lower=c(0,0,0), upper=c(1,1,1))
	# The worst direction detected is nearly proportional to (2,-1,2)
	# (type "?OA131" for explanations about this linear orthogonal array)
print(rss$worst.dir)

# Now, scramble this design
# X <- (OA131 + matrix(runif(49*3, 49, 3)))/7
# or load the design obtained this way
data(OA131_scrambled)
OA2 <- OA131_scrambled

# no feature is detected by the 2D RSS:
rss <- rss2d(OA2, lower=c(0,0,0), upper=c(1,1,1))    
# 4 clusters are detected by the 3D RSS:
rss <- rss3d(OA2, lower=c(0,0,0), upper=c(1,1,1))	

	
	# Defect detection of 8D Sobol sequences
	# All triplets of dimensions are tried to detect the worst defect 
	# (according to the specified goodness-of-fit statistic).
	# requires randtoolbox library to generate the Sobol sequence
## Not run: 
##D library(randtoolbox)
##D d <- 8
##D n <- 10*d
##D rss <- rss3d(design=sobol(n=n, dim=d), lower=rep(0,d), upper=rep(1,d))
## End(Not run)



cleanEx()
nameEx("runif.faure")
### * runif.faure

flush(stderr()); flush(stdout())

### Name: runif.faure
### Title: Low discrepancy sequence : Faure
### Aliases: runif.faure
### Keywords: design

### ** Examples

f <- runif.faure(20,2)
plot(f$design,xlim=c(0,1),ylim=c(0,1))



cleanEx()
nameEx("straussDesign")
### * straussDesign

flush(stderr()); flush(stdout())

### Name: straussDesign
### Title: Designs based on Strauss process
### Aliases: straussDesign
### Keywords: design

### ** Examples

# Strauss-Gibbs designs in dimension 2 (n=20 points)
S1 <- straussDesign(n=20,dimension=2,RND=0.2)
plot(S1$design,xlim=c(0,1),ylim=c(0,1))
theta <- seq(0,2*pi,by =2*pi/(100 - 1))
for(i in 1:S1$n){
   lines(S1$design[i,1]+S1$radius/2*cos(theta),
	   S1$design[i,2]+S1$radius/2*sin(theta),col='red')
}
# 2D-Strauss design
S2 <- straussDesign(n=20,dimension=2,RND=0.2,NMC=200,
	constraints1D=0,alpha=0,repulsion=0.01)

plot(S2$design,xlim=c(0,1),ylim=c(0,1))

# 2D-Strauss designs with constraints on the axis
S3 <- straussDesign(n=20,dimension=2,RND=0.18,NMC=200,
	constraints1D=1,alpha=0.5,repulsion=0.1,repulsion1D=0.01)

plot(S3$design,xlim=c(0,1),ylim=c(0,1))
rug(S3$design[,1],side=1)
rug(S3$design[,2],side=2)



cleanEx()
nameEx("wspDesign")
### * wspDesign

flush(stderr()); flush(stdout())

### Name: wspDesign
### Title: WSP algorithm
### Aliases: wspDesign
### Keywords: design

### ** Examples

dimension <- 2
n <- 100
X <- matrix(runif(n*dimension),n,dimension)
m=wspDesign(X,0.1)
plot(m)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
