CPP_smooth.manifold.FEM.basis<-function(locations, observations, FEMbasis, covariates = NULL, incidence_matrix = NULL, ndim, mydim, BC = NULL, opt_method, lambdas = NULL, initial_lambda = NULL, nrealizations = 100)
{
  # C++ function for manifold works with vectors not with matrices
  
  FEMbasis$mesh$triangles=c(t(FEMbasis$mesh$triangles))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))
  
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
 
  FEMbasis$mesh$triangles=FEMbasis$mesh$triangles-1
  
  
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }
  
  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  { 
    BC$BC_indices<-as.vector(BC$BC_indices)-1
    
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }
  
  if(is.null(lambdas))
  {
    lambdas<-vector(length=0)
  }else
  {
    lambdas<-as.vector(lambdas)
  }
  
  if(is.null(initial_lambda))
  {
    initial_lambda = 0
  }
  
  opt_method<-as.vector(opt_method)
  
  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  data <- as.vector(observations)
  storage.mode(observations) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nnodes) <- "integer"
  storage.mode(FEMbasis$mesh$ntriangles) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values)  <- "double"
  storage.mode(opt_method) <- "integer"
  storage.mode(lambdas) <- "double"
  storage.mode(initial_lambda) <- "double"
  storage.mode(nrealizations) <- "integer"
  
  ## Call C++ function
  bigsol <- .Call("regression_Laplace", locations, data, FEMbasis$mesh, FEMbasis$mesh$order, mydim, ndim, covariates,
                  incidence_matrix, BC$BC_indices, BC$BC_values, opt_method, lambdas, initial_lambda, nrealizations, PACKAGE = "fdaPDE")
  
  return(bigsol)
}

CPP_eval.manifold.FEM = function(FEM, locations, incidence_matrix, redundancy, ndim, mydim)
{
  FEMbasis = FEM$FEMbasis
  
  # C++ function for manifold works with vectors not with matrices
  
  FEMbasis$mesh$triangles=c(t(FEMbasis$mesh$triangles))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))
  
  #NEW: riporto shift indici in R
  FEMbasis$mesh$triangles=FEMbasis$mesh$triangles-1
  
  # Imposing types, this is necessary for correct reading from C++
  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  coeff <- as.matrix(FEM$coeff)
  storage.mode(coeff) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(locations) <- "double"
  storage.mode(redundancy) <- "integer"
  
  #Calling the C++ function "eval_FEM_fd" in RPDE_interface.cpp
  evalmat = matrix(0,max(nrow(locations),nrow(incidence_matrix)),ncol(coeff))
  for (i in 1:ncol(coeff)){
    evalmat[,i] <- .Call("eval_FEM_fd", FEMbasis$mesh, locations, incidence_matrix, coeff[,i],
                         FEMbasis$order, redundancy, mydim, ndim, PACKAGE = "fdaPDE")
  }
  
  #Returning the evaluation matrix
  evalmat
}

