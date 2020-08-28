##########################################
############## TEST SCRIPTS FOR THE SIMULATIONS IN THE REPORT ###############
##########################################

#if not intsalled, install the following libraries
#install.packages("tictoc")
#install.packages("writexl")

library(fdaPDE)
library(tictoc)
library(writexl)
      
#build logarithmic sequence function
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

####### 2D ########

#### Test 1: square domain ####
#            locations = nodes 
#            laplacian    
#            no covariates
#            no BC
#            order FE = 1


N=25  #Number of repetitions
nums=seq(1,N,by=1)
data1 <- data.frame(tests=nums)     #Create a dataframe for the values to compute


x = seq(0,1, length.out = 41)
y = x
locations = expand.grid(x,y)
  
  mesh = create.mesh.2D(locations)  #refine mesh for other cases
  #refinements of mesh used: 
  #no refine 1681 nodes
  # mesh=refine.mesh.2D(mesh, maximum_area = 2e-4)   #3281 nodes
  # mesh=refine.mesh.2D(mesh, maximum_area = 5e-5)   #15741 nodes
  # mesh=refine.mesh.2D(mesh, maximum_area = 2.5e-5)   #30971 nodes

  plot(mesh)

nodes=mesh$nodes
nnodes=dim(mesh$nodes)[1]

FEMbasis=create.FEM.basis(mesh)

# Test function
f = function(x, y, z = 1){
  coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
  sin(2*pi*(coe(y,1)*x*cos(z-2)-y*sin(z-2)))*cos(2*pi*(coe(y,1)*x*cos(z-2+pi/2)+coe(x,1)*y*sin((z-2)*pi/2)))
}

# Exact solution (pointwise at nodes)
sol_exact=f(mesh$nodes[,1], mesh$nodes[,2])

image(FEM(sol_exact, FEMbasis))
#dev.off()


#rgl.snapshot(filename="bst32.png",fmt="png")   #save images 
  

# Set smoothing parameter
lambda= lseq(10^-7,10^2,length.out = 20)




####Test1.1: grid of 20 values GCV exact

#Set important 
MSE_nocov=NULL
mse=NULL
time_tot=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)

for (i in 1:N){
graphics.off()

# Add error to simulate data
set.seed(10000*i)
ran=range(sol_exact)
data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
                       
image(output_CPP1$fit.FEM)



colnames(GCV_param1) = c("edf", "sd", "GCV")
#compute time for each call to smooth.FEM
tic()
  output_CPP1<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda, optimization = "grid", loss_function = "GCV", DOF_evaluation = "exact")
t=toc()
time=t$toc-t$tic
time_tot=c(time_tot,time)
  #plot(output_CPP1$fit.FEM)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSE=sum((sol_exact-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE_nocov=c(MSE_nocov, MSE)
    ##STORE GCV PARAMETERS
    GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
    GCV_param1[i,2]=output_CPP1$solution$estimated_sd
    GCV_param1[i,3]=output_CPP1$optimization$GCV
}


RMSE_nocov=sqrt(MSE_nocov)
data1$grid_exact_RMSE=RMSE_nocov  #save RMSE
data1$grid_exact_time=time_tot #Save time for each computation
data1$grid_exact_lambdas=selected_lambda #save best lambda
data1$grid_exact_GCV=GCV_param1[,3]  #save best GCV value
data1$grid_exact_dof=GCV_param1[,1]  #save dof of best GCV
data1$grid_exact_sd=GCV_param1[,2]   #save sd



# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test1.2: Newton exact with exact GCV


MSE_nocov=NULL
mse=NULL
time_tot=NULL
selected_lambda=rep(0,N)
GCV_param1 = matrix(nrow=N, ncol=3)

for (i in 1:N){
  graphics.off()
  
  # Add error to simulate data
  set.seed(10000*i)
  ran=range(sol_exact)
  data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  colnames(GCV_param1) = c("edf", "sd", "GCV")

  #compute time
  tic()
  output_CPP1<-smooth.FEM(observations=data, FEMbasis=FEMbasis, optimization = "newton", loss_function = "GCV", DOF_evaluation = "exact", stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  #plot(output_CPP1$fit.FEM)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSE=sum((sol_exact-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE_nocov=c(MSE_nocov, MSE)
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
}


RMSE_nocov=sqrt(MSE_nocov) 
data1$newton_exact_RMSE=RMSE_nocov #save RMSE
data1$newton_exact_time=time_tot  #Save time for each computation
data1$newton_exact_mse=mse
data1$newton_exact_lambdas=selected_lambda  #save best lambda
data1$newton_exact_GCV=GCV_param1[,3]    #save best GCV value
data1$newton_exact_dof=GCV_param1[,1]    #save dof of best GCV
data1$newton_exact_sd=GCV_param1[,2]       #save sd



# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

####Test1.3: Newton finite differences with exact GCV

MSE_nocov=NULL
mse=NULL
time_tot=NULL
selected_lambda=rep(0,N)  
GCV_param1 = matrix(nrow=N, ncol=3)
for (i in 1:N){
  graphics.off()
  
  # Add error to simulate data
  set.seed(10000*i)
  ran=range(sol_exact)
  data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  #image(output_CPP$fit.FEM)
  
  

  colnames(GCV_param1) = c("edf", "sd", "GCV")
 
  #compute time
  tic()
  output_CPP1<-smooth.FEM(observations=data, FEMbasis=FEMbasis, optimization = "newton_fd", loss_function = "GCV", DOF_evaluation = "exact", stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)

  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSE=sum((sol_exact-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE_nocov=c(MSE_nocov, MSE)
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
}


RMSE_nocov=sqrt(MSE_nocov)
data1$newton_fd_exact_RMSE=RMSE_nocov #save RMSE
data1$newton_fd_exact_time=time_tot   #save time for each computation
data1$newton_fd_exact_lambdas=selected_lambda  #save best lambda
data1$newton_fd_exact_GCV=GCV_param1[,3]     #save best GCV value
data1$newton_fd_exact_dof= GCV_param1[,1]     #save dof of best GCV
data1$newton_fd_exact_sd=GCV_param1[,2]       #save sd



# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
  boxplot(time_tot)
  
####Test1.4: Newton finite differences with stochastic GCV

  
MSE_nocov=NULL
mse=NULL
time_tot=NULL
selected_lambda=rep(0,N)
GCV_param1 = matrix(nrow=N, ncol=3)
for (i in 1:N){
  graphics.off()
  
  # Add error to simulate data
  set.seed(10000*i)
  ran=range(sol_exact)
  data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  

  colnames(GCV_param1) = c("edf", "sd", "GCV")

  #compute time 
  tic()
  output_CPP1<-smooth.FEM(observations=data, FEMbasis=FEMbasis, optimization = "newton_fd", loss_function = "GCV", DOF_evaluation = "stochastic", seed=10000*i, stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)

  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSE=sum((sol_exact-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE_nocov=c(MSE_nocov, MSE)
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
}
  

RMSE_nocov=sqrt(MSE_nocov)
data1$newton_fd_stoch_RMSE=RMSE_nocov #save RMSE
data1$newton_fd_stoch_time=time_tot   #save time for each computation

data1$newton_fd_stoch_lambdas=selected_lambda    #save best lambda
data1$newton_fd_stoch_GCV=GCV_param1[,3]         #save best GCV value
data1$newton_fd_stoch_dof=GCV_param1[,1]         #save dof of best GCV
data1$newton_fd_stoch_sd=GCV_param1[,2]          #save sd
 

# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

####Test1.5: grid of 20 values GCV stochastic


MSE_nocov=NULL
mse=NULL
time_tot=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
  
for (i in 1:N){
  graphics.off()
  
  # Add error to simulate data
  set.seed(10000*i)
  ran=range(sol_exact)
  data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  #image(output_CPP$fit.FEM)
  
  
  
  colnames(GCV_param1) = c("edf", "sd", "GCV")
  
  #compute time
  tic()
  output_CPP1<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda, optimization = "grid", loss_function = "GCV", DOF_evaluation = "stochastic", seed=10000*i)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  #plot(output_CPP1$fit.FEM)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSE=sum((sol_exact-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE_nocov=c(MSE_nocov, MSE)
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}

RMSE_nocov=sqrt(MSE_nocov)
data1$grid_stoch_RMSE=RMSE_nocov  #save RNSE
data1$grid_stoch_time=time_tot    #save time for each computation
data1$grid_stoch_lambdas=selected_lambda   #save best lambda
data1$grid_stoch_GCV=GCV_param1[,3]          #save best GCV value
data1$grid_stoch_dof=GCV_param1[,1]          #save dof of best GCV 
data1$grid_stoch_sd=GCV_param1[,2]              #save sd


# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

# Mean RMSE
MSE_hat=mean(MSE_nocov)
mse_hat=mean(mse)
time_hat=mean(time_tot)



library("writexl")  #save dataframe in the current directory
write_xlsx(data1,"data_refine_rett_tol0_005.xlsx")



library(fdaPDE)
library(tictoc)




#build logarithmic sequence
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}


#### Test 2: c-shaped domain ####
#            locations != nodes
#            laplacian
#            with covariates
#            no BC
#            order FE = 1

N=25  #number of repetitions
nums=seq(1,N,by=1)   
data2 <- data.frame(tests=nums)   #crate dataframe


graphics.off()

data(horseshoe2D)

mesh=create.mesh.2D(nodes=horseshoe2D$boundary_nodes, 
                    segments = horseshoe2D$boundary_segments)
locations=refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)

plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")


FEMbasis=create.FEM.basis(mesh)

ndata = nrow(locations)
nnodes=dim(mesh$nodes)[1]

nodes=mesh$nodes
f_ex=fs.test(nodes[,1],nodes[,2])   #exact f 


image(FEM(f_ex, FEMbasis))

####Test2.1: grid of 20 values GCV exact


# Set smoothing parameter
lambda= 10^seq(-7,2,length.out=20)

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL


GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  
  # Create covariates
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(locations[,1])
  W=cbind(cov1, cov2)
  beta_exact=c(2,-1)
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2]) + 2*cov1 - cov2
  
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  #compute time
  tic()
  output_CPP1<-smooth.FEM(locations = locations, observations=data, 
                          covariates = cbind(cov1, cov2), 
                          FEMbasis=FEMbasis, lambda=lambda, optimization = "grid", loss_function = "GCV", DOF_evaluation = "exact")
  t=toc()
  
  
 # image(output_CPP1$fit.FEM)
 # rgl.snapshot(filename="newtst_0.005.png",fmt="png")  #save images

  time=t$toc-t$tic
  time_tot=c(time_tot,time)

  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=sum((f_ex-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE=sum((as.vector(DatiEsatti)-(output_CPP1$solution$z_hat))^2)/dim(locations)[1]
  MSE_nocov=c(MSE_nocov, MSEnp)  #MSE used
  MSE_g=c(MSE_g,MSE)
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
}


RMSE_nocov=sqrt(MSE_nocov)
RMSE_g=sqrt(MSE_g)
data2$grid_exact_RMSE=RMSE_nocov  #RMSE used in the report
data2$grid_exact_time=time_tot  #time
data2$grid_exact_lambdas=selected_lambda  #best lambda
data2$grid_exact_GCV=GCV_param1[,3]    #best GCV
data2$grid_exact_dof=GCV_param1[,1]    #dof of best GCV
data2$grid_exact_sd=GCV_param1[,2]     #sd of best GCV
data2$grid_exact_RMSEg=RMSE_g          #RMSE with covariates, not used
data2$grid_exact_beta1=beta1           #beta hat 1
data2$grid_exact_beta2=beta2           #beta hat 2

# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test2.2: Newton exact with exact GCV

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")

beta1=NULL
beta2=NULL

for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  
  # Create covariates
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(locations[,1])
  W=cbind(cov1, cov2)
  beta_exact=c(2,-1)
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2]) + 2*cov1 - cov2
  
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  #time
  tic()
  output_CPP1<-smooth.FEM(locations = locations, observations=data, 
                          covariates = cbind(cov1, cov2),
                          FEMbasis=FEMbasis, optimization = "newton", loss_function = "GCV", DOF_evaluation = "exact", stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=sum((f_ex-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE=sum((as.vector(DatiEsatti)-(output_CPP1$solution$z_hat))^2)/dim(locations)[1]
  MSE_nocov=c(MSE_nocov, MSEnp)
  MSE_g=c(MSE_g,MSE)
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}


RMSE_nocov=sqrt(MSE_nocov)   #RMSE used 
RMSE_g=sqrt(MSE_g)
data2$newton_exact_RMSE=RMSE_nocov   #RMSE used
data2$newton_exact_time=time_tot    #time
data2$newton_exact_lambdas=selected_lambda  #best lambda
data2$newton_exact_GCV=GCV_param1[,3]   #best GCV
data2$newton_exact_dof=GCV_param1[,1]   #dof
data2$newton_exact_sd=GCV_param1[,2]    #sd
data2$newton_exact_RMSEg=RMSE_g         #RMSE with covriates, not used in the test   
data2$newton_exact_beta1=beta1         #beta hat 1
data2$newton_exact_beta2=beta2         #beta hat 2

# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

####Test2.3: Newton finite differences with exact GCV
MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  
  # Create covariates
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(locations[,1])
  W=cbind(cov1, cov2)
  beta_exact=c(2,-1)
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2]) + 2*cov1 - cov2

  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  tic()
  output_CPP1<-smooth.FEM(locations = locations, observations=data, 
                          covariates = cbind(cov1, cov2),
                          FEMbasis=FEMbasis, optimization = "newton_fd", loss_function = "GCV", DOF_evaluation = "exact", stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=sum((f_ex-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE=sum((as.vector(DatiEsatti)-(output_CPP1$solution$z_hat))^2)/dim(locations)[1]
  mse=c(mse, output_CPP1$solution$rmse)
  MSE_nocov=c(MSE_nocov, MSEnp)
  MSE_g=c(MSE_g,MSE)
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
}


RMSE_nocov=sqrt(MSE_nocov)
RMSE_g=sqrt(MSE_g)
data2$newton_fd_exact_RMSE=RMSE_nocov
data2$newton_fd_exact_time=time_tot
data2$newton_fd_exact_mse=mse
data2$newton_fd_exact_lambdas=selected_lambda
data2$newton_fd_exact_GCV=GCV_param1[,3]
data2$newton_fd_exact_dof=GCV_param1[,1]
data2$newton_fd_exact_sd=GCV_param1[,2]
data2$newton_fd_exact_RMSEg=RMSE_g
data2$newton_fd_exact_beta1=beta1
data2$newton_fd_exact_beta2=beta2



# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)
x11()
boxplot(mse)
x11()
boxplot(RMSE_g)
# Mean RMSE
MSE_hat=mean(MSE_nocov)
mse_hat=mean(mse)
time_hat=mean(time_tot)



####Test1.4: grid of 15 values stochastic

# Set smoothing parameter
lambda= 10^seq(-7,2,length.out=20)

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  
  # Create covariates
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(locations[,1])
  W=cbind(cov1, cov2)
  beta_exact=c(2,-1)
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2]) + 2*cov1 - cov2

  
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  output_CPP1<-smooth.FEM(locations = locations, observations=data, 
                          covariates = cbind(cov1, cov2),
                          FEMbasis=FEMbasis, lambda=lambda, optimization = "grid", loss_function = "GCV", DOF_evaluation = "stochastic", seed=10000*i)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=sum((f_ex-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE=sum((as.vector(DatiEsatti)-(output_CPP1$solution$z_hat))^2)/dim(locations)[1]
  MSE_nocov=c(MSE_nocov, MSEnp)
  MSE_g=c(MSE_g,MSE)
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  
}

RMSE_nocov=sqrt(MSE_nocov)
RMSE_g=sqrt(MSE_g)
data2$grid_stoch_RMSE=RMSE_nocov  #RMSE used
data2$grid_stoch_time=time_tot    #total time
data2$grid_stoch_lambdas=selected_lambda    #optimal lambda
data2$grid_stoch_GCV=GCV_param1[,3]        #GCV optimal
data2$grid_stoch_dof=GCV_param1[,1]        # dof
data2$grid_stoch_sd=GCV_param1[,2]         #sd
data2$grid_stoch_RMSEg=RMSE_g             #not used
data2$grid_stoch_beta1=beta1              #beta hat 1
data2$grid_stoch_beta2=beta2              #beta hat 2 

# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test1.5: Newton finite differences with stochastic GCV
# Set smoothing parameter

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  
  # Create covariates
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(locations[,1])
  W=cbind(cov1, cov2)
  beta_exact=c(2,-1)
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2]) + 2*cov1 - cov2

  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  output_CPP1<-smooth.FEM(locations = locations, observations=data, 
                          covariates = cbind(cov1, cov2), 
                          FEMbasis=FEMbasis, optimization = "newton_fd", loss_function = "GCV", DOF_evaluation = "stochastic", seed=10000*i, stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=sum((f_ex-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE=sum((as.vector(DatiEsatti)-(output_CPP1$solution$z_hat))^2)/dim(locations)[1]
  MSE_nocov=c(MSE_nocov, MSEnp)
  MSE_g=c(MSE_g,MSE)
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  
}



RMSE_nocov=sqrt(MSE_nocov)  
RMSE_g=sqrt(MSE_g)
data2$newton_fd_stoch_RMSE=RMSE_nocov   #RMSE used
data2$newton_fd_stoch_time=time_tot     #time
data2$newton_fd_stoch_lambdas=selected_lambda  #best lambda
data2$newton_fd_stoch_GCV=GCV_param1[,3]       #best GCV
data2$newton_fd_stoch_dof=GCV_param1[,1]      #dof
data2$newton_fd_stoch_sd=GCV_param1[,2]       #sd
data2$newton_fd_stoch_RMSEg=RMSE_g            #not used
data2$newton_fd_stoch_beta1=beta1             #beta hat 1
data2$newton_fd_stoch_beta2=beta2             #beta hat 2


# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

graphics.off()

library("writexl")
write_xlsx(data2,"data_correct_c_noref_tol0_05.xlsx")   #save dataframe in the current directory





#build logarithmic sequence function
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

#### Test 3: c-shaped domain ####
#            locations != nodes
#            laplacian
#            without covariates
#            no BC
#            order FE = 1

N=25  #number of repetitions
nums=seq(1,N,by=1)   
data3 <- data.frame(tests=nums)   #crate dataframe


graphics.off()

data(horseshoe2D)

mesh=create.mesh.2D(nodes=horseshoe2D$boundary_nodes, 
                    segments = horseshoe2D$boundary_segments)
locations=refine.mesh.2D(mesh, maximum_area = 0.05)$nodes
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)

plot(mesh, pch = ".")
points(locations, pch = 16, col = "red")


FEMbasis=create.FEM.basis(mesh)

ndata = nrow(locations)
nnodes=dim(mesh$nodes)[1]

nodes=mesh$nodes
f_ex=fs.test(nodes[,1],nodes[,2])   #exact f 


image(FEM(f_ex, FEMbasis))

####Test3.1: grid of 20 values GCV exact


# Set smoothing parameter
lambda= 10^seq(-7,2,length.out=20)

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL


GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  
  # Create covariates
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(locations[,1])
  W=cbind(cov1, cov2)
  beta_exact=c(2,-1)
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2]) 
  
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  #compute time
  tic()
  output_CPP1<-smooth.FEM(locations = locations, observations=data, 
                        
                          FEMbasis=FEMbasis, lambda=lambda, optimization = "grid", loss_function = "GCV", DOF_evaluation = "exact")
  t=toc()
  
  
  #image(output_CPP1$fit.FEM)
  # rgl.snapshot(filename="newtst_0.005.png",fmt="png")  #save images
  
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=sum((f_ex-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE=sum((as.vector(DatiEsatti)-(output_CPP1$solution$z_hat))^2)/dim(locations)[1]
  MSE_nocov=c(MSE_nocov, MSEnp)  #MSE used
  MSE_g=c(MSE_g,MSE)
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
}


RMSE_nocov=sqrt(MSE_nocov)
RMSE_g=sqrt(MSE_g)
data3$grid_exact_RMSE=RMSE_nocov  #RMSE used in the report
data3$grid_exact_time=time_tot  #time
data3$grid_exact_lambdas=selected_lambda  #best lambda
data3$grid_exact_GCV=GCV_param1[,3]    #best GCV
data3$grid_exact_dof=GCV_param1[,1]    #dof of best GCV
data3$grid_exact_sd=GCV_param1[,2]     #sd of best GCV
data3$grid_exact_RMSEg=RMSE_g          #RMSE with covariates, not used
data3$grid_exact_beta1=beta1           #beta hat 1
data3$grid_exact_beta2=beta2           #beta hat 2

# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test3.2: Newton exact with exact GCV

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")

beta1=NULL
beta2=NULL

for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  
  # Create covariates
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(locations[,1])
  W=cbind(cov1, cov2)
  beta_exact=c(2,-1)
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2])
  
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  #time
  tic()
  output_CPP1<-smooth.FEM(locations = locations, observations=data, 
                     
                          FEMbasis=FEMbasis, optimization = "newton", loss_function = "GCV", DOF_evaluation = "exact", stop_criterion_tol = 0.005)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=sum((f_ex-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE=sum((as.vector(DatiEsatti)-(output_CPP1$solution$z_hat))^2)/dim(locations)[1]
  MSE_nocov=c(MSE_nocov, MSEnp)
  MSE_g=c(MSE_g,MSE)
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}


RMSE_nocov=sqrt(MSE_nocov)   #RMSE used 
RMSE_g=sqrt(MSE_g)
data3$newton_exact_RMSE=RMSE_nocov   #RMSE used
data3$newton_exact_time=time_tot    #time
data3$newton_exact_lambdas=selected_lambda  #best lambda
data3$newton_exact_GCV=GCV_param1[,3]   #best GCV
data3$newton_exact_dof=GCV_param1[,1]   #dof
data3$newton_exact_sd=GCV_param1[,2]    #sd
data3$newton_exact_RMSEg=RMSE_g         #RMSE with covriates, not used in the test   
data3$newton_exact_beta1=beta1         #beta hat 1
data3$newton_exact_beta2=beta2         #beta hat 2

# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

####Test3.3: Newton finite differences with exact GCV
MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  
  # Create covariates
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(locations[,1])
  W=cbind(cov1, cov2)
  beta_exact=c(2,-1)
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2]) 
  
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  tic()
  output_CPP1<-smooth.FEM(locations = locations, observations=data, 
                         
                          FEMbasis=FEMbasis, optimization = "newton_fd", loss_function = "GCV", DOF_evaluation = "exact", stop_criterion_tol = 0.005)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  #plot(output_CPP1$fit.FEM)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=sum((f_ex-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE=sum((as.vector(DatiEsatti)-(output_CPP1$solution$z_hat))^2)/dim(locations)[1]
  mse=c(mse, output_CPP1$solution$rmse)
  MSE_nocov=c(MSE_nocov, MSEnp)
  MSE_g=c(MSE_g,MSE)
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
}


RMSE_nocov=sqrt(MSE_nocov)
RMSE_g=sqrt(MSE_g)
data3$newton_fd_exact_RMSE=RMSE_nocov
data3$newton_fd_exact_time=time_tot
data3$newton_fd_exact_mse=mse
data3$newton_fd_exact_lambdas=selected_lambda
data3$newton_fd_exact_GCV=GCV_param1[,3]
data3$newton_fd_exact_dof=GCV_param1[,1]
data3$newton_fd_exact_sd=GCV_param1[,2]
data3$newton_fd_exact_RMSEg=RMSE_g
data3$newton_fd_exact_beta1=beta1
data3$newton_fd_exact_beta2=beta2



# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)
x11()
boxplot(mse)
x11()
boxplot(RMSE_g)
# Mean RMSE
MSE_hat=mean(MSE_nocov)
mse_hat=mean(mse)
time_hat=mean(time_tot)



####Test3.4: grid of 20 values stochastic

# Set smoothing parameter
lambda= 10^seq(-7,2,length.out=20)

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  
  # Create covariates
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(locations[,1])
  W=cbind(cov1, cov2)
  beta_exact=c(2,-1)
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2]) 
  
  for(ind in 1:100){
    points(locations[which(round((DatiEsatti-min(DatiEsatti))/(max(DatiEsatti)-min(DatiEsatti))*100)==ind),1],locations[which(round((DatiEsatti-min(DatiEsatti))/(max(DatiEsatti)-min(DatiEsatti))*100)==ind),2],
           col=heat.colors(100)[ind], pch=16)
  }
  
  
  
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  output_CPP1<-smooth.FEM(locations = locations, observations=data, 
                         
                          FEMbasis=FEMbasis, lambda=lambda, optimization = "grid", loss_function = "GCV", DOF_evaluation = "stochastic", seed=10000*i)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=sum((f_ex-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE=sum((as.vector(DatiEsatti)-(output_CPP1$solution$z_hat))^2)/dim(locations)[1]
  MSE_nocov=c(MSE_nocov, MSEnp)
  MSE_g=c(MSE_g,MSE)
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  
}

RMSE_nocov=sqrt(MSE_nocov)
RMSE_g=sqrt(MSE_g)
data3$grid_stoch_RMSE=RMSE_nocov  #RMSE used
data3$grid_stoch_time=time_tot    #total time
data3$grid_stoch_lambdas=selected_lambda    #optimal lambda
data3$grid_stoch_GCV=GCV_param1[,3]        #GCV optimal
data3$grid_stoch_dof=GCV_param1[,1]        # dof
data3$grid_stoch_sd=GCV_param1[,2]         #sd
data3$grid_stoch_RMSEg=RMSE_g             #not used
data3$grid_stoch_beta1=beta1              #beta hat 1
data3$grid_stoch_beta2=beta2              #beta hat 2 

# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test3.5: Newton finite differences with stochastic GCV
# Set smoothing parameter

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  
  # Create covariates
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(locations[,1])
  W=cbind(cov1, cov2)
  beta_exact=c(2,-1)
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test(locations[,1], locations[,2]) 
  
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  output_CPP1<-smooth.FEM(locations = locations, observations=data, 
                    
                          FEMbasis=FEMbasis, optimization = "newton_fd", loss_function = "GCV", DOF_evaluation = "stochastic", seed=10000*i, stop_criterion_tol = 0.005)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=sum((f_ex-output_CPP1$fit.FEM$coeff)^2)/nnodes
  MSE=sum((as.vector(DatiEsatti)-(output_CPP1$solution$z_hat))^2)/dim(locations)[1]
  MSE_nocov=c(MSE_nocov, MSEnp)
  MSE_g=c(MSE_g,MSE)
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  
}



RMSE_nocov=sqrt(MSE_nocov)  
RMSE_g=sqrt(MSE_g)
data3$newton_fd_stoch_RMSE=RMSE_nocov   #RMSE used
data3$newton_fd_stoch_time=time_tot     #time
data3$newton_fd_stoch_lambdas=selected_lambda  #best lambda
data3$newton_fd_stoch_GCV=GCV_param1[,3]       #best GCV
data3$newton_fd_stoch_dof=GCV_param1[,1]      #dof
data3$newton_fd_stoch_sd=GCV_param1[,2]       #sd
data3$newton_fd_stoch_RMSEg=RMSE_g            #not used
data3$newton_fd_stoch_beta1=beta1             #beta hat 1
data3$newton_fd_stoch_beta2=beta2             #beta hat 2


# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

graphics.off()

library("writexl")
write_xlsx(data3,"data_correct_c_noref_tol0_005no_cov.xlsx")   #save dataframe in the current directory





rm(list=ls())
graphics.off()
#build logarithmic sequence
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}


#### Test 4: square domain ####
#            locations in nodes
#            PDE
#            no covariates
#            no BC
#            order FE = 2
####### 2D ########

N=25
nums=seq(1,N,by=1)     #cretae datafrane
data4 <- data.frame(tests=nums)

x = seq(0,1, length.out = 11)
y = x
locations = expand.grid(x,y)

mesh = create.mesh.2D(locations, order = 2)
plot(mesh)
FEMbasis=create.FEM.basis(mesh)

# Test function
a1=1
a2=4
z<-function(p){  
  a1*sin(2*pi*p[,1])*cos(2*pi*p[,2])+a2*sin(3*pi*p[,1])}

# Exact solution (pointwise at nodes)
sol_exact=z(mesh$nodes)
image(FEM(sol_exact, FEMbasis))

DatiEsatti=z(locations)
ndati = length(DatiEsatti)

# Set smoothing parameter
lambda= 10^seq(-7,2,length.out = 20)

# Set PDE parameters
PDE_parameters = list(K = matrix(c(1,0,0,4), nrow = 2), b = c(0,0), c = 0)

####Test4.1: grid of 20 values GCV exact

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL


GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")

for (i in 1:N){
  
  # Add error to simulate data
  #set.seed(7893475)
  set.seed(10000*i)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(ndati, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  tic()
  
  output_CPP1<-smooth.FEM(observations=data, 
                          FEMbasis=FEMbasis, 
                          lambda=lambda,
                          PDE_parameters=PDE_parameters,
                          optimization = "grid", loss_function = "GCV", DOF_evaluation = "exact")
  t=toc()
  
  image(output_CPP1$fit.FEM)
  #rgl.snapshot(filename="nex0.05.png",fmt="png")     #take image of the field
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  plot(output_CPP1$fit.FEM)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)
  MSE_nocov=c(MSE_nocov, MSEnp)
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}


RMSE_nocov=sqrt(MSE_nocov)
data4$grid_exact_RMSE=RMSE_nocov   #RMSE
data4$grid_exact_time=time_tot     #Time
data4$grid_exact_lambdas=selected_lambda    #optimal lambda
data4$grid_exact_GCV=GCV_param1[,3]      #GCV optimal
data4$grid_exact_dof=GCV_param1[,2]     #dof of optimal GCV
data4$grid_exact_sd=GCV_param1[,1]      #sd standard deviation of optimal GCV

# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test4.2: Newton exact with exact GCV

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")

beta1=NULL
beta2=NULL

for (i in 1:N){
  
  # Add error to simulate data
  #set.seed(7893475)
  set.seed(10000*i)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(ndati, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  
  output_CPP1<-smooth.FEM(observations=data, 
                          FEMbasis=FEMbasis, 
                          PDE_parameters=PDE_parameters,
                          optimization = "newton", loss_function = "GCV", DOF_evaluation = "exact", stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)
  MSE_nocov=c(MSE_nocov, MSEnp)

    ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
  
  
}


RMSE_nocov=sqrt(MSE_nocov)

data4$newton_exact_RMSE=RMSE_nocov    #RMSE
data4$newton_exact_time=time_tot      #Time
data4$newton_exact_lambdas=selected_lambda    #lambda opt
data4$newton_exact_GCV=GCV_param1[,3]        #GCV opt
data4$newton_exact_dof=GCV_param1[,1]        #dof of GCV opt
data4$newton_exact_sd=GCV_param1[,2]        #sd of GCV opt

# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test4.3: Newton finite differences with exact GCV
# Set smoothing parameter

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  
  # Add error to simulate data
  #set.seed(7893475)
  set.seed(10000*i)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(ndati, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  
  output_CPP1<-smooth.FEM(observations=data, 
                          FEMbasis=FEMbasis, 
                          PDE_parameters=PDE_parameters,
                          optimization = "newton_fd", loss_function = "GCV", DOF_evaluation = "exact", stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)
  MSE_nocov=c(MSE_nocov, MSEnp)
  
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}


RMSE_nocov=sqrt(MSE_nocov)     
data4$newton_fd_exact_RMSE=RMSE_nocov   #RMSE
data4$newton_fd_exact_time=time_tot     #Time 
data4$newton_fd_exact_lambdas=selected_lambda    #best lambda
data4$newton_fd_exact_GCV=GCV_param1[,3]    #GCV opt
data4$newton_fd_exact_dof=GCV_param1[,1]     #dof of GCV opt
data4$newton_fd_exact_sd=GCV_param1[,2]      #sd of GCV opt

# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test4.4: grid of 20 values GCV stochastic

# Set smoothing parameter
lambda= 10^seq(-7,2,length.out=20)

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  
  # Add error to simulate data
  #set.seed(7893475)
  set.seed(10000*i)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(ndati, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  
  output_CPP1<-smooth.FEM(observations=data, 
                          FEMbasis=FEMbasis, 
                          lambda=lambda,
                          PDE_parameters=PDE_parameters,
                          optimization = "grid", loss_function = "GCV", DOF_evaluation = "stochastic", seed = 10000*i)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)
  MSE_nocov=c(MSE_nocov, MSEnp)
 
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}

RMSE_nocov=sqrt(MSE_nocov)

data4$grid_stoch_RMSE=RMSE_nocov   #RMSE
data4$grid_stoch_time=time_tot     #Time
data4$grid_stoch_lambdas=selected_lambda   #opt lambda
data4$grid_stoch_GCV=GCV_param1[,3]    #opt GCV
data4$grid_stoch_dof=GCV_param1[,1]     #diof of opt GCV
data4$grid_stoch_sd=GCV_param1[,2]     #sd of opt GCV


# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test4.5: Newton finite differences with stochastic GCV

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  
  # Add error to simulate data
  #set.seed(7893475)
  set.seed(10000*i)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(ndati, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  
  output_CPP1<-smooth.FEM(observations=data, 
                          FEMbasis=FEMbasis, 
                          lambda=lambda,
                          PDE_parameters=PDE_parameters,
                          optimization = "newton_fd", loss_function = "GCV", DOF_evaluation = "stochastic", seed=10000*i, stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)
  MSE_nocov=c(MSE_nocov, MSEnp)

  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
  
}


RMSE_nocov=sqrt(MSE_nocov)                  
data4$newton_fd_stoch_RMSE=RMSE_nocov    #RMSE
data4$newton_fd_stoch_time=time_tot      #Time
data4$newton_fd_stoch_lambdas=selected_lambda    #opt lambda
data4$newton_fd_stoch_GCV=GCV_param1[,3]       #opt GCV
data4$newton_fd_stoch_dof=GCV_param1[,1]      #dof of opt GCV
data4$newton_fd_stoch_sd=GCV_param1[,2]     #sd of opt GCV



# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

graphics.off()

library("writexl")
write_xlsx(data4,"data_pde_0_05tol.xlsx")                #save dataframe in the current directory






#build logarithmic sequence
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

#####The next test is used to shot the monotonicity of GCV 

####### 2D ########

N=25
nums=seq(1,N,by=1)
data5 <- data.frame(tests=nums)   #create dataframe



#### Test 5: quasicircular domain ####
#            areal observations
#            PDE space varying
#            no covariates
#            with BC
#            order FE = 1



data(quasicircle2Dareal)
mesh = quasicircle2Dareal$mesh
incidence_matrix = quasicircle2Dareal$incidence_matrix
DatiEsatti = quasicircle2Dareal$data

plot(mesh)
FEMbasis = create.FEM.basis(mesh)

image(FEM(DatiEsatti, FEMbasis))

#rgl.snapshot(filename="bst32.png",fmt="png")  #take snapshot of the fieds

# Set smoothing parameter
lambda =  10^seq(-7,2,length.out=20)

# Set BC
BC = NULL
BC$BC_indices = which(mesh$nodesmarkers == 1)
BC$BC_values = rep(0,length(BC$BC_indices))

# Set sv-PDE parameters
R = 2.8
K1 = 0.1
K2 = 0.2
beta = 0.5

K_func<-function(points)
{
  output = array(0, c(2, 2, nrow(points)))
  for (i in 1:nrow(points))
    output[,,i] = 10*rbind(c(points[i,2]^2 + K1*points[i,1]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2),
                             (K1-1)*points[i,1]*points[i,2]),
                           c((K1-1)*points[i,1]*points[i,2],
                             points[i,1]^2 + K1*points[i,2]^2 + K2*(R^2 - points[i,1]^2 - points[i,2]^2)))
  output
}

b_func<-function(points)
{
  output = array(0, c(2, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = 10*beta*c(points[i,1],points[i,2])
  output
}

c_func<-function(points)
{
  rep(c(0), nrow(points))
}


#### Forcing term != 0
# forcing function != 0
u_func<-function(points)
{
  output = array(0, c(1, nrow(points)))
  for (i in 1:nrow(points))
    output[,i] = -ifelse((points[i,1]^2+points[i,2]^2)<1,100,0)
  output
}

# plot the forcing funcion
xgrid=seq(from=-3,to=3,by=0.1)
ygrid=seq(from=-3,to=3,by=0.1)
image(xgrid,ygrid,matrix(u_func(expand.grid(xgrid,ygrid)),nrow=length(xgrid),ncol=length(ygrid),byrow=FALSE),main='forcing function',asp=1)
lines(mesh$nodes[1:50,])

PDE_parameters = list(K = K_func, b = b_func, c = c_func, u = u_func)
nloc=7 #7 exact data

####Test5.1: grid of 20 values GCV exact

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL


GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  # Add error to simulate data
  #set.seed(5839745)
  data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.05*(max(DatiEsatti)-min(DatiEsatti)))
  
  
  data = data + 5
  
  # Set new value for the BC
  BC$BC_values = rep(5,length(BC$BC_indices))

  tic()
  
  
  output_CPP1<-smooth.FEM(observations=data, 
                          incidence_matrix = incidence_matrix,
                          FEMbasis=FEMbasis, 
                          lambda=lambda,
                          BC = BC, 
                          PDE_parameters = PDE_parameters, 
                          optimization='grid', DOF_evaluation='exact', loss_function='GCV')
  t=toc()
  
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  image(output_CPP1$fit.FEM)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((DatiEsatti-output_CPP1$solution$z_hat)^2)
  MSE_nocov=c(MSE_nocov, MSEnp)

  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}


RMSE_nocov=sqrt(MSE_nocov)

data5$grid_exact_RMSE=RMSE_nocov               #RMSE
data5$grid_exact_time=time_tot                 #Time
data5$grid_exact_lambdas=selected_lambda       #lambda opt
data5$grid_exact_GCV=GCV_param1[,3]            #GCV opt
data5$grid_exact_dof=GCV_param1[,2]           #dof of GCV opt
data5$grid_exact_sd=GCV_param1[,1]           #sd of GCV opt


# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test5.2: Newton finite differences with exact GCV
# Set smoothing parameter

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  # Add error to simulate data
  #set.seed(5839745)
  data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.05*(max(DatiEsatti)-min(DatiEsatti)))
  
  
  data = data + 5
  
  # Set new value for the BC
  BC$BC_values = rep(5,length(BC$BC_indices))
  
  
  tic()
  
  
  output_CPP1<-smooth.FEM(observations=data, 
                          incidence_matrix = incidence_matrix,
                          FEMbasis=FEMbasis, 
                          
                          BC = BC, 
                          PDE_parameters = PDE_parameters, 
                          optimization='newton_fd', DOF_evaluation='exact', loss_function='GCV', stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((DatiEsatti-output_CPP1$solution$z_hat)^2)

  MSE_nocov=c(MSE_nocov, MSEnp)
 
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}


RMSE_nocov=sqrt(MSE_nocov)

data5$newton_fd_exact_RMSE=RMSE_nocov  #RMSE
data5$newton_fd_exact_time=time_tot    #Time
data5$newton_fd_exact_lambdas=selected_lambda   #lambda opt
data5$newton_fd_exact_GCV=GCV_param1[,3]        #GCV opt
data5$newton_fd_exact_dof=GCV_param1[,1]         #df of opt GCV
data5$newton_fd_exact_sd=GCV_param1[,2]         #sd of opt GCV



# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)



####Test5.3: grid of 20 values GCV stochastic

# Set smoothing parameter
lambda= 10^seq(-7,2,length.out=20)

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  # Add error to simulate data
  #set.seed(5839745)
  data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.05*(max(DatiEsatti)-min(DatiEsatti)))
  
  
  data = data + 5
  
  # Set new value for the BC
  BC$BC_values = rep(5,length(BC$BC_indices))
  
  
  tic()
  
  
  output_CPP1<-smooth.FEM(observations=data, 
                          incidence_matrix = incidence_matrix,
                          FEMbasis=FEMbasis, 
                          lambda=lambda,
                          BC = BC, 
                          PDE_parameters = PDE_parameters, 
                          optimization='grid', DOF_evaluation='stochastic', loss_function='GCV', seed=10000*i)
  t=toc()
  
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((DatiEsatti-output_CPP1$solution$z_hat)^2)
  MSE_nocov=c(MSE_nocov, MSEnp)
 
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}

RMSE_nocov=sqrt(MSE_nocov)
data5$grid_stoch_RMSE=RMSE_nocov          #RMSE
data5$grid_stoch_time=time_tot            #Time
data5$grid_stoch_lambdas=selected_lambda  #Lambda opt
data5$grid_stoch_GCV=GCV_param1[,3]      #GCV of lambda opt
data5$grid_stoch_dof=GCV_param1[,1]      #dof of lambda opt
data5$grid_stoch_sd=GCV_param1[,2]       #sd of lambda opt

# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test5.4: Newton finite differences with stochastic GCV

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  set.seed(10000*i)
  # Add error to simulate data
  #set.seed(5839745)
  data = DatiEsatti + rnorm(length(DatiEsatti), sd = 0.05*(max(DatiEsatti)-min(DatiEsatti)))
  
  
  data = data + 5
  
  # Set new value for the BC
  BC$BC_values = rep(5,length(BC$BC_indices))
  
  
  tic()
  
  
  output_CPP1<-smooth.FEM(observations=data, 
                          incidence_matrix = incidence_matrix,
                          FEMbasis=FEMbasis, 
                          
                          BC = BC, 
                          PDE_parameters = PDE_parameters, 
                          optimization='newton_fd', DOF_evaluation='stochastic', loss_function='GCV', seed=10000*i, stop_criterion_tol = 0.05)
  t=toc()
  
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((DatiEsatti-output_CPP1$solution$z_hat)^2)
  MSE_nocov=c(MSE_nocov, MSEnp)

  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
  
}



RMSE_nocov=sqrt(MSE_nocov)

data5$newton_fd_stoch_RMSE=RMSE_nocov           #RMSE
data5$newton_fd_stoch_time=time_tot             #Time
data5$newton_fd_stoch_lambdas=selected_lambda   #lambda opt
data5$newton_fd_stoch_GCV=GCV_param1[,3]        #GCV opt
data5$newton_fd_stoch_dof=GCV_param1[,1]        #dof of GCV opt
data5$newton_fd_stoch_sd=GCV_param1[,2]         #sd of GCV opt

# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

library("writexl")
write_xlsx(data5,"data_areal.xlsx")                 #save dataframe in the current directory





#build logarithmic sequence function
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

#### Test 6: sphere domain ####
#            locations != nodes 
#            with covariates
#            no BC
#            order FE = 1
# ONLY beta is used (only one covariate)
N=25
nums=seq(1,N,by=1)
data6 <- data.frame(tests=nums)       #create dataframe

graphics.off()

data(sphere2.5D)
mesh = sphere2.5D

FEMbasis=create.FEM.basis(mesh)

# Test function 
f = function(x, y, z){
  phi = (asin(y/sqrt(x^2+y^2)))
  theta = acos(z/sqrt(x^2+y^2+z^2))
  # rho = 1
  
  sin(4*(1/2*sin(theta)*exp(-sin(theta)^2)+1)*theta)*cos(2*(1/2*cos(phi)*exp(-cos(phi)^2)+1)*phi)
}

# Exact solution (pointwise at nodes)
sol_exact=f(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3])
plot(FEM(sol_exact, FEMbasis))

#rgl.snapshot(filename="bst.png",fmt="png")                snapshot of the field


# Set smoothing parameter
lambda = 10^seq(-7,2,length.out = 20)

####Test 6.1: grid of 20 values GCV exact

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL
  

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  set.seed(10000*i)
  # Add error to simulate data
  
  # Generate data locations on the sphere
  #set.seed(598944)
  ndata = 500
  locations = matrix(rnorm(ndata*3), ncol = 3)
  locations = locations/sqrt(locations[,1]^2 + locations[,2]^2 + locations[,3]^2)
  
  # Generate covariate and data
  cov1 = runif(ndata, min = -1, max = 1)
  DatiEsatti = f(locations[,1],locations[,2], locations[,3]) + cov1
  # set.seed(7893475)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(ndata, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  # Project locations on the mesh
  projected_locations = projection.points.2.5D(mesh, locations)

  tic()
  
  output_CPP1<-smooth.FEM(observations=data, locations = projected_locations,
                          covariates = cov1,
                          FEMbasis=FEMbasis, lambda=lambda,
                          optimization='grid', DOF_evaluation='exact', loss_function='GCV')
  t=toc()
  plot(output_CPP1$fit.FEM)

  time=t$toc-t$tic
  time_tot=c(time_tot,time)

  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)
  
  MSE_nocov=c(MSE_nocov, MSEnp)

  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}



RMSE_nocov=sqrt(MSE_nocov)

data6$grid_exact_RMSE=RMSE_nocov  #RMSE
data6$grid_exact_time=time_tot   #time

data6$grid_exact_lambdas=selected_lambda    #optimal lambda
data6$grid_exact_GCV=GCV_param1[,3]         #opt GCV
data6$grid_exact_dof=GCV_param1[,2]         #dof of opt GCV
data6$grid_exact_sd=GCV_param1[,1]         #sd of opt GCV

data6$grid_exact_beta1=beta1             #beta hat 1
data6$grid_exact_beta2=beta2             #beta hat 2

# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)





####Test6.2: Newton exact with exact GCV

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")

beta1=NULL
beta2=NULL

for (i in 1:N){
  # Add error to simulate data
  # set.seed(7893475)
  set.seed(10000*i)
  
  # Generate data locations on the sphere
  #set.seed(598944)
  ndata = 500
  locations = matrix(rnorm(ndata*3), ncol = 3)
  locations = locations/sqrt(locations[,1]^2 + locations[,2]^2 + locations[,3]^2)
  
  # Generate covariate and data
  cov1 = runif(ndata, min = -1, max = 1)
  DatiEsatti = f(locations[,1],locations[,2], locations[,3]) + cov1
  
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(ndata, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  # Project locations on the mesh
  projected_locations = projection.points.2.5D(mesh, locations)
  
  
  tic()
  
  output_CPP1<-smooth.FEM(observations=data, locations = projected_locations,
                          covariates = cov1,
                          FEMbasis=FEMbasis,
                          optimization='newton', DOF_evaluation='exact', loss_function='GCV', stop_criterion_tol = 0.05)
  t=toc()

  time=t$toc-t$tic
  time_tot=c(time_tot,time)
 
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)

 
  MSE_nocov=c(MSE_nocov, MSEnp)
 
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
  
  
}


RMSE_nocov=sqrt(MSE_nocov)                    #RMSE
data6$newton_exact_time=time_tot              #time
data6$newton_exact_lambdas=selected_lambda    #lambda OPT
data6$newton_exact_GCV=GCV_param1[,3]          #GCV opt
data6$newton_exact_dof=GCV_param1[,1]          #dof of GCV opt
data6$newton_exact_sd=GCV_param1[,2]           #sd of GCV opt
 
data6$newton_exact_beta1=beta1                 #beta hat 1 opt
data6$newton_exact_beta2=beta2

# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test6.3: Newton finite differences with exact GCV

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  # set.seed(7893475)
  
  set.seed(10000*i)
  # Generate data locations on the sphere
  #set.seed(598944)
  ndata = 500
  locations = matrix(rnorm(ndata*3), ncol = 3)
  locations = locations/sqrt(locations[,1]^2 + locations[,2]^2 + locations[,3]^2)
  
  # Generate covariate and data
  cov1 = runif(ndata, min = -1, max = 1)
  DatiEsatti = f(locations[,1],locations[,2], locations[,3]) + cov1
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(ndata, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  # Project locations on the mesh
  projected_locations = projection.points.2.5D(mesh, locations)
  
  
  tic()
  
  output_CPP1<-smooth.FEM(observations=data, locations = projected_locations,
                          covariates = cov1,
                          FEMbasis=FEMbasis,
                          optimization='newton_fd', DOF_evaluation='exact', loss_function='GCV', stop_criterion_tol = 0.05)
  t=toc()
 
  time=t$toc-t$tic
  time_tot=c(time_tot,time)

  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)
 
  
  MSE_nocov=c(MSE_nocov, MSEnp)

  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}


RMSE_nocov=sqrt(MSE_nocov)

data6$newton_fd_exact_RMSE=RMSE_nocov          #RMSE
data6$newton_fd_exact_time=time_tot           #time

data6$newton_fd_exact_lambdas=selected_lambda  #lambda opt
data6$newton_fd_exact_GCV=GCV_param1[,3]       #GCV opt
data6$newton_fd_exact_dof=GCV_param1[,1]        #dof of opt GCV
data6$newton_fd_exact_sd=GCV_param1[,2]         #sd of opt GCV

data6$newton_fd_exact_beta1=beta1               #beta hat 1
data6$newton_fd_exact_beta2=beta2



# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test6.4: grid of 20 values GCV stochastic

# Set smoothing parameter
lambda= 10^seq(-8,3,length.out=20)

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  # set.seed(7893475)
  set.seed(10000*i)
  
  
  # Generate data locations on the sphere
  #set.seed(598944)
  ndata = 500
  locations = matrix(rnorm(ndata*3), ncol = 3)
  locations = locations/sqrt(locations[,1]^2 + locations[,2]^2 + locations[,3]^2)
  
  # Generate covariate and data
  cov1 = runif(ndata, min = -1, max = 1)
  DatiEsatti = f(locations[,1],locations[,2], locations[,3]) + cov1
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(ndata, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  # Project locations on the mesh
  projected_locations = projection.points.2.5D(mesh, locations)
  
  
  tic()
  
  output_CPP1<-smooth.FEM(observations=data, locations = projected_locations,
                          covariates = cov1,
                          FEMbasis=FEMbasis, lambda=lambda,
                          optimization='grid', DOF_evaluation='stochastic', loss_function='GCV',  seed=10000*i)
  t=toc()
  
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  #plot(output_CPP1$fit.FEM)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)

 
  MSE_nocov=c(MSE_nocov, MSEnp)

  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
  
}

RMSE_nocov=sqrt(MSE_nocov)

data6$grid_stoch_RMSE=RMSE_nocov        #RMSE
data6$grid_stoch_time=time_tot          #time

data6$grid_stoch_lambdas=selected_lambda    #lambda opt
data6$grid_stoch_GCV=GCV_param1[,3]         #GCV opt
data6$grid_stoch_dof=GCV_param1[,1]        #dof of GCV opt
data6$grid_stoch_sd=GCV_param1[,2]         #sd of GCV opt

data6$grid_stoch_beta1=beta1               #beta hat 1
data6$grid_stoch_beta2=beta2

# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test6.5: Newton finite differences with stochastic GCV
# Set smoothing parameter

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Add error to simulate data
  # set.seed(7893475)
  set.seed(10000*i)
  
  # Generate data locations on the sphere
  #set.seed(598944)
  ndata = 500
  locations = matrix(rnorm(ndata*3), ncol = 3)
  locations = locations/sqrt(locations[,1]^2 + locations[,2]^2 + locations[,3]^2)
  
  # Generate covariate and data
  cov1 = runif(ndata, min = -1, max = 1)
  DatiEsatti = f(locations[,1],locations[,2], locations[,3]) + cov1
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(ndata, mean=0, sd=0.05*abs(ran[2]-ran[1]))
  
  # Project locations on the mesh
  projected_locations = projection.points.2.5D(mesh, locations)
  
  
  tic()
  
  output_CPP1<-smooth.FEM(observations=data, locations = projected_locations,
                          covariates = cov1,
                          FEMbasis=FEMbasis,
                          optimization='newton_fd', DOF_evaluation='stochastic', loss_function='GCV', seed=10000*i, stop_criterion_tol = 0.05)
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)
  MSE_nocov=c(MSE_nocov, MSEnp)

  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
  
}



RMSE_nocov=sqrt(MSE_nocov)

data6$newton_fd_stoch_RMSE=RMSE_nocov             #RMSE
data6$newton_fd_stoch_time=time_tot               #time

data6$newton_fd_stoch_lambdas=selected_lambda     #lambda opt
data6$newton_fd_stoch_GCV=GCV_param1[,3]          #GCV opt
data6$newton_fd_stoch_dof=GCV_param1[,1]          #dof of GCV opt
data6$newton_fd_stoch_sd=GCV_param1[,2]           #sd of GCV opt

# Boxplot of RMSE
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

graphics.off()

library("writexl")
write_xlsx(data6,"data_sphere_0_05tol_solo_stoch.xlsx")                 #save dataframe in the current directory






#### Test 7: c-shaped domain ####
#            locations = nodes
#            with covariates
#            no BC
#            order FE = 1



rm(list=ls())
graphics.off()
#build logarithmic sequence
lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}



N=25
nums=seq(1,N,by=1)
data7 <- data.frame(tests=nums)    #create dataframe


graphics.off()

data(horseshoe3D)

mesh = horseshoe3D

FEMbasis=create.FEM.basis(mesh)

ndata = FEMbasis$nbasis
sol_exact=fs.test.3D(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]) 

#rgl.snapshot(filename="3d/exact.png",fmt="png")             #snapshot of the fields

plot(output_CPP1$fit.FEM)


####Test1.1: grid of 15 values exact


MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL


GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")

#set smoothing parameters 
lambda = 10^seq(-7,3,length.out = 20)

for (i in 1:N){
  # Create covariates
  #set.seed(509875)
  set.seed(10000*i)
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(mesh$nodes[,1])
  
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test.3D(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]) + 2*cov1 - cov2
  #plot(FEM(DatiEsatti, FEMbasis))
  
  # Add error to simulate data
  #set.seed(543663)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean   =0, sd=0.05*abs(ran[2]-ran[1]))

  tic()
  output_CPP1<-smooth.FEM(observations=data, 
                          covariates = cbind(cov1, cov2),
                          FEMbasis=FEMbasis, lambda=lambda,
                          optimization='grid', DOF_evaluation='exact', loss_function='GCV')
  t=toc()
  time=t$toc-t$tic
  time_tot=c(time_tot,time)

  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)

  MSE_nocov=c(MSE_nocov, MSEnp)

  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
}



RMSE_nocov=sqrt(MSE_nocov)

data7$grid_exact_RMSE=RMSE_nocov        #RMSE
data7$grid_exact_time=time_tot           #time
 
data7$grid_exact_lambdas=selected_lambda  #lambda opt
data7$grid_exact_GCV=GCV_param1[,3]         #GCV opt
data7$grid_exact_dof=GCV_param1[,2]       #dof of opt GCV
data7$grid_exact_sd=GCV_param1[,1]        #sd of opt GCV

data7$grid_exact_beta1=beta1             #beta 1 hat 
data7$grid_exact_beta2=beta2             #beta 2 hat

# Boxplots 
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test7.2: Newton exact with exact GCV

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")

beta1=NULL
beta2=NULL

for (i in 1:N){
  # Create covariates
  #set.seed(509875)
  set.seed(10000*i)
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(mesh$nodes[,1])
  
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test.3D(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]) + 2*cov1 - cov2
  #plot(FEM(DatiEsatti, FEMbasis))
  
  # Add error to simulate data
  #set.seed(543663)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean   =0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  output_CPP1<-smooth.FEM(observations=data, 
                          covariates = cbind(cov1, cov2),
                          FEMbasis=FEMbasis,
                          optimization='newton', DOF_evaluation='exact', loss_function='GCV', stop_criterion_tol = 0.05)
  t=toc()
 
  time=t$toc-t$tic
  time_tot=c(time_tot,time)
  
  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)


  MSE_nocov=c(MSE_nocov, MSEnp)
  
  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
  
  
  
}


RMSE_nocov=sqrt(MSE_nocov)

data7$newton_exact_RMSE=RMSE_nocov   #RMSE
data7$newton_exact_time=time_tot     #time

data7$newton_exact_lambdas=selected_lambda #lambda opt
data7$newton_exact_GCV=GCV_param1[,3]      #GCV opt
data7$newton_exact_dof=GCV_param1[,1]      #dof of GCV opt
data7$newton_exact_sd=GCV_param1[,2]       #sd of GCV opt

data7$newton_exact_beta1=beta1            #beta1 hat
data7$newton_exact_beta2=beta2            #beta2 hat

# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

####Test7.3: Newton finite differences with exact GCV


MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Create covariates
  #set.seed(509875)
  set.seed(10000*i)
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(mesh$nodes[,1])
  
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test.3D(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]) + 2*cov1 - cov2
  #plot(FEM(DatiEsatti, FEMbasis))
  
  # Add error to simulate data
  #set.seed(543663)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean   =0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  output_CPP1<-smooth.FEM(observations=data, 
                          covariates = cbind(cov1, cov2),
                          FEMbasis=FEMbasis,
                          optimization='newton_fd', DOF_evaluation='exact', loss_function='GCV', stop_criterion_tol = 0.05)
  t=toc()
  
  time=t$toc-t$tic
  time_tot=c(time_tot,time)

  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)

  MSE_nocov=c(MSE_nocov, MSEnp)

  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
}


RMSE_nocov=sqrt(MSE_nocov)

data7$newton_fd_exact_RMSE=RMSE_nocov        #RNSE
data7$newton_fd_exact_time=time_tot         #time
data7$newton_fd_exact_lambdas=selected_lambda   #lambda opt
data7$newton_fd_exact_GCV=GCV_param1[,3]        #GCV opt
data7$newton_fd_exact_dof=GCV_param1[,1]         #dof of GCV opt
data7$newton_fd_exact_sd=GCV_param1[,2]         #sd of GCV opt

data7$newton_fd_exact_beta1=beta1                #beta1 hat
data7$newton_fd_exact_beta2=beta2                #beta2 hat



# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test7.4: grid of 20 values GCV stochastic

# Set smoothing parameter
lambda= 10^seq(-7,3,length.out=20)

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Create covariates
  #set.seed(509875)
  set.seed(10000*i)
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(mesh$nodes[,1])
  
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test.3D(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]) + 2*cov1 - cov2
  #plot(FEM(DatiEsatti, FEMbasis))
  
  # Add error to simulate data
  #set.seed(543663)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean   =0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  output_CPP1<-smooth.FEM(observations=data, 
                          covariates = cbind(cov1, cov2),
                          FEMbasis=FEMbasis, lambda=lambda,
                          optimization='grid', DOF_evaluation='stochastic', loss_function='GCV', seed=10000*i)
  t=toc()
 
  time=t$toc-t$tic
  time_tot=c(time_tot,time)

  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)
  
  
  MSE_nocov=c(MSE_nocov, MSEnp)

  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof[lambda_index]
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
  
}

RMSE_nocov=sqrt(MSE_nocov)

data7$grid_stoch_RMSE=RMSE_nocov       #RMSE
data7$grid_stoch_time=time_tot         #time
 
data7$grid_stoch_lambdas=selected_lambda  #lambda_opt
data7$grid_stoch_GCV=GCV_param1[,3]       #GCV opt
data7$grid_stoch_dof=GCV_param1[,1]       #dof of GCV opt
data7$grid_stoch_sd=GCV_param1[,2]        #sd of GCV opt

data7$grid_stoch_beta1=beta1              #beta1 hat
data7$grid_stoch_beta2=beta2              #beta2 hat

# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)


####Test7.5: Newton finite differences with stochastic GCV
# Set smoothing parameter

MSE_nocov=NULL
MSE_g=NULL
mse=NULL
time_tot=NULL
beta1=NULL
beta2=NULL

GCV_param1 = matrix(nrow=N, ncol=3)
selected_lambda=rep(0,N)
colnames(GCV_param1) = c("edf", "sd", "GCV")


for (i in 1:N){
  # Create covariates
  #set.seed(509875)
  set.seed(10000*i)
  cov1=rnorm(ndata, mean = 1, sd = 2)
  cov2=sin(mesh$nodes[,1])
  
  # Exact solution (pointwise at nodes)
  DatiEsatti=fs.test.3D(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]) + 2*cov1 - cov2
  #plot(FEM(DatiEsatti, FEMbasis))
  
  # Add error to simulate data
  #set.seed(543663)
  ran=range(DatiEsatti)
  data = DatiEsatti + rnorm(length(DatiEsatti), mean   =0, sd=0.05*abs(ran[2]-ran[1]))
  
  
  tic()
  output_CPP1<-smooth.FEM(observations=data, 
                          covariates = cbind(cov1, cov2),
                          FEMbasis=FEMbasis,
                          optimization='newton_fd', DOF_evaluation='stochastic', loss_function='GCV', stop_criterion_tol = 0.005, seed=10000*i)
  t=toc()
  
  time=t$toc-t$tic
  time_tot=c(time_tot,time)

  selected_lambda[i]=output_CPP1$optimization$lambda_solution
  lambda_index=output_CPP1$optimization$lambda_position
  MSEnp=mean((sol_exact-output_CPP1$fit.FEM$coeff)^2)

  MSE_nocov=c(MSE_nocov, MSEnp)

  beta1=c(beta1, output_CPP1$solution$beta[1])
  beta2=c(beta2,output_CPP1$solution$beta[2])
  ##STORE GCV PARAMETERS
  GCV_param1[i,1]=output_CPP1$optimization$dof
  GCV_param1[i,2]=output_CPP1$solution$estimated_sd
  GCV_param1[i,3]=output_CPP1$optimization$GCV
}



RMSE_nocov=sqrt(MSE_nocov)
data7$newton_fd_stoch_RMSE=RMSE_nocov   #RMSE
data7$newton_fd_stoch_time=time_tot     #Time

data7$newton_fd_stoch_lambdas=selected_lambda        #lambda opt
data7$newton_fd_stoch_GCV=GCV_param1[,3]            #GCV opt
data7$newton_fd_stoch_dof=GCV_param1[,1]           #dof of GCV opt
data7$newton_fd_stoch_sd=GCV_param1[,2]             #sd of GCV opt
 
data7$newton_fd_stoch_beta1=beta1                   #beta 1 hat
data7$newton_fd_stoch_beta2=beta2                   #beta 2 hat 


# Boxplots
x11()
boxplot(RMSE_nocov)
x11()
boxplot(time_tot)

graphics.off()

#Close all the windows open
while (rgl.cur() > 0) { rgl.close() }

library("writexl")
write_xlsx(data7,"data_3d.xlsx")  #save dataframe in the current directory







