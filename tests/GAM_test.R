
library(fdaPDE)

  #### square 2D (basic case)
  # setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace/fdaPDE_ultimate/data")
  rm(list=ls())
  graphics.off()
  
  x = seq(0,1, length.out = 5)
  y = x
  locations = expand.grid(x,y)
  
  mesh = create.mesh.2D(locations)
  plot(mesh)
  
  nnodes=dim(mesh$nodes)[1]
  
  FEMbasis=create.FEM.basis(mesh)
  


    set.seed(5847947)
    # set.seed(42)
    
    a1=runif(1,min=-1.5,max=1.5)
    a2=runif(1,min=-1.5,max=1.5)
    
    z<-function(p)
    {
      
      a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1]) - 2
      
    }
    
    # Exact solution (pointwise at nodes)
    sol_exact=rep(0,dim(mesh$nodes)[1])
    for(i in 1: dim(mesh$nodes)[1])
      sol_exact[i]=z(mesh$nodes[i,])
    
    ran=range(sol_exact) 
    
  
  
  # Set smoothing parameter ---------------------------
  
    
    link<-function(x){-1/x}
    inv.link<-function(x){-1/x}
    
    count = 0;
    desmat=matrix(0,nrow=nnodes,ncol=2)
    mu=numeric(length=nnodes)
    scale.param = 1
    
    while(any(mu<=0))
    {
      count <- count+1
      # The seed is set in such a way that if count <10, which is very likely, then
      # the seeds are different for all the simulations (which is highly desirable!)
      set.seed(42 + count)
      desmat[,1]=rbeta(nnodes,shape1=1.5,shape2=2)
      desmat[,2]=rbeta(nnodes,shape1=3,shape2=2)+1
      beta1=-2/5
      beta2=3/10
      betas_truth = c(beta1,beta2)
      param=sol_exact+beta1*desmat[,1]+beta2*desmat[,2]
      mu<-inv.link(param)
      #if(all(mu>=0)) seeds.check[sim] <- seeds[sim] + count
    }
    
    response <- rgamma(nnodes, shape=mu/scale.param, scale=scale.param)
  
  #to give the dofs for  lambda= c(10^-4,10^-3)
  dofs=  matrix(c(19.10020, 13.97923), nrow = 2)
    
    
  GCVFLAG=T
  GCVMETHODFLAG='Exact'
  lambda= c(10^-4,10^-3)
  
  #mu_guessed <- rep(1,nnodes) + response
  
    
  output_CPP <- smooth.FEM(observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                                    lambda = lambda, max.steps=10, family="gamma", mu0=NULL, scale.param=scale.param,  DOF_matrix = dofs, loss_function = "GCV")
  