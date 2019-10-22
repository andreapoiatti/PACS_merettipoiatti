###########################
## test smooth.FEM.basis ##
###########################

library(fdaPDE)
source('pointsFEM.R')

GCVFLAG=TRUE  # do not compute GCV (default)

# RECOMMENDATION: after testing without GCV computation, repeat test with both GCV methods

# GCVMETHODFLAG=1 for exact GCV
# GCVMETHODFLAG=2 for stochastic GCV (default)

# load the C-shaped mesh
load("meshC.RData")
FEMbasis = create.FEM.basis(mesh)

# exact solution f: Ramsay's test function

fs.test <- function (x, y, r0 = 0.1, r = 0.5, l = 3, b = 1, exclude = TRUE) 
{
  
  q <- pi * r/2
  a <- d <- x * 0
  
  ind <- x >= 0 & y > 0
  a[ind] <- q + x[ind]
  d[ind] <- y[ind] - r
  
  ind <- x >= 0 & y <= 0
  a[ind] <- (-q - x[ind])
  d[ind] <- -r - y[ind]
  
  ind <- x < 0
  a[ind] <- -atan(y[ind]/x[ind]) * r 
  d[ind] <-( sqrt(x[ind]^2 + y[ind]^2) - r )* (y[ind]/r0*(as.numeric(abs(y[ind])<=r0 & x[ind]>-0.5))+(as.numeric(abs(y[ind])>r0 || x[ind]<(-0.5))))
  
  ind <- abs(d) > r - r0 | (x > l & (x - l)^2 + d^2 > (r - r0)^2)
  
  f <- a * b + d^2
  
  if (exclude) 
    f[ind] <- NA
  
  attr(f, "exclude") <- ind
  f
}

### locations at mesh nodes

# generate exact data from the test function
dati_esatti = fs.test(mesh$nodes[,1], mesh$nodes[,2], exclude = FALSE)

# plot of exact solution
image(FEM(dati_esatti,FEMbasis = FEMbasis))
x11()
points.2D.data(SpacePoints = mesh$nodes, Data = dati_esatti)

# add gaussian noise to data

sd = 5*abs(max(dati_esatti) - min(dati_esatti))/100

set.seed(5847947) # MEMO: do not forget the seed to replicate the simulation!

data = dati_esatti + rnorm(length(dati_esatti), sd = sd)
x11()
points.2D.data(SpacePoints = mesh$nodes, Data = data)

# choose value of lambda parameter
lambda = 10^-2

# solution
solution = smooth.FEM.basis(observations = data, FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG)
plot(solution$fit.FEM)

points=eval.FEM(solution$fit.FEM, locations = mesh$nodes)
write.table(points,file="smoothFEMbasis2D_nod_nocov.txt")


### locations different from nodes

set.seed(5847947) 

loc1=cbind(runif(50,min=0.5,max=2.5),runif(50,min=0.1,max=0.9)) #braccio superiore
loc2=cbind(runif(50,min=0.5,max=2.5),runif(50,min=-0.9,max=-0.1)) # braccio inferiore
loc3=cbind(runif(50,min=-0.7,max=-0.1),runif(50,min=-0.5,max=0.5)) #circonferenza grande
noditest=mesh$nodes[1:50,]# alcune oss coincidenti con nodi

oss=rbind(loc1,loc2,loc3,noditest)

dati_esatti2 = fs.test(oss[,1], oss[,2], exclude = FALSE)

x11()
points.2D.data(SpacePoints = oss, Data = dati_esatti2)


sd2 = 5*abs(max(dati_esatti2) - min(dati_esatti2))/100

data2 = dati_esatti2 + rnorm(length(dati_esatti2), sd = sd2)
x11()
points.2D.data(SpacePoints = oss, Data = data2)


lambda = 0.01


solution = smooth.FEM.basis(locations=oss,observations = data2, FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG)
plot(solution$fit.FEM)

points=eval.FEM(solution$fit.FEM,locations=mesh$nodes)
write.table(points,"smoothFEMbasis2D_nonod_nocov.txt")

### locations at nodes, with covariates

cov=cbind(rnorm(dim(mesh$nodes)[1],mean=0.5,sd=0.01), rnorm(dim(mesh$nodes)[1], mean=2.5,sd=0.05))
solution = smooth.FEM.basis(observations = data, covariates=cov,FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG, GCVmethod = 1)
plot(solution$fit.FEM)
solution$beta
#[1,] 
#[2,]  

x11()
points.2D.data(SpacePoints = mesh$nodes,Data=cov[,1])
x11()
points.2D.data(SpacePoints = mesh$nodes,Data=cov[,2])
points=eval.FEM(solution$fit.FEM,locations=mesh$nodes)
write.table(points,"smoothFEMbasis2D_nod_cov.txt")
#write.table(solution$beta)

### locations different from nodes, with covariates

cov=cbind(rnorm(dim(oss)[1],mean=0.5,sd=0.01), rnorm(dim(oss)[1], mean=2.5,sd=0.05))
solution = smooth.FEM.basis(locations=oss,observations = data2, covariates=cov,FEMbasis = FEMbasis, lambda = lambda, GCV=GCVFLAG, GCVmethod = 1)
plot(solution$fit.FEM)
x11()
points.2D.data(SpacePoints = oss,Data=cov[,1])
x11()
points.2D.data(SpacePoints = oss,Data=cov[,2])
solution$beta
# [1,]  
# [2,] 

x11()
points.2D.data(SpacePoints = oss,Data=cov[,1])
x11()
points.2D.data(SpacePoints = oss,Data=cov[,2])
points=eval.FEM(solution$fit.FEM,locations=oss)
write.table(points,"smoothFEMbasis2D_nonod_cov.txt")
