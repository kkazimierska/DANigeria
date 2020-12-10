# Input 
# spdf - spatial polygon or spatial polygon data frame, that we create the polygon 

# max.edge - vector of the max.edge for the initial mesh based on the boundary
# and for the mesh based on the initial mesh locations

# cutoff - vector of the cutoff for the initial mesh based on the boundary
# and for the mesh based on the initial mesh locations

# offset - offset for the mesh based on the initial mesh locations

# Output 
# mesh traingulation 

create_mesh <- function(spdf, max.edge, cutoff, offset, plot = T){
  
  ng.bound = unionSpatialPolygons(
    spdf, IDs = rep(1, length(spdf))) 
  
  ### define an initial mesh inner domain
  mesh0 <- inla.mesh.2d(
    boundary = ng.bound, 
    max.edge = max.edge[1], 
    cutoff = cutoff[1])   
  mesh0$n #Parana was 5 200 nodes
  
  
  mesh <- inla.mesh.2d(
    mesh0$loc[, 1:2], 
    max.edge = max.edge[2], 
    offset = offset,   
    cutoff = cutoff[2])    

  if(plot == TRUE){
    plot(mesh, asp=1, lwd=2)
    plot(spdf, add=TRUE, border=gray(.7))
  } else {
    
  }
  
  return(mesh)
  
}



