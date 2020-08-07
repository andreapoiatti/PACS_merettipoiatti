library(fdaPDE)
data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations
time_locations = seq(0,1,length.out = 5)

mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
plot(mesh)
#mesh=refine.mesh.2D(mesh, maximum_area = 1e-3)
#plot(mesh)
space_time_locations = cbind(rep(time_locations,each=nrow(mesh$nodes)),
                             rep(mesh$nodes[,1],5),rep(mesh$nodes[,2],5))

FEMbasis = create.FEM.basis(mesh)
lambdaS = c(10^-1)
lambdaT = c(10^-1)
set.seed(1234)
data = fs.test(space_time_locations[,2],
               space_time_locations[,3])*cos(pi*space_time_locations[,1]) +
  rnorm(nrow(space_time_locations), sd = 0.5)
data = matrix(data, nrow = nrow(mesh$nodes), ncol = length(time_locations), byrow = TRUE)

solution = smooth.FEM.time(observations = data, time_locations = time_locations,
                           FEMbasis = FEMbasis, lambdaS = lambdaS, lambdaT = lambdaT, GCV=TRUE, FLAG_PARABOLIC = FALSE, GCVmethod = "Exact")
plot(solution$fit.FEM)
