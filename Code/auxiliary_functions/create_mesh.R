# Input 
# spdf - spatial polygon or spatial polygon data frame, that we create the polygon 

# Output 
# mesh traingulation 

create_mesh <- function(spdf, plot = T ){
  
  ng.bound = unionSpatialPolygons(
    spdf, IDs = rep(1, nrow(spdf))) 
  
  ext <- extent(coordinates(spdf))
  range.xy = c(ext@xmax-ext@xmin, ext@ymax-ext@ymin)
  
  
  ### define an initial mesh inner domain
  mesh0 <- inla.mesh.2d(
    boundary=ng.bound, 
    max.edge=floor(0.078*mean(range.xy)/4), 
    cutoff=floor(0.078*mean(range.xy)/4)/2)   
  mesh0$n #Parana was 5 200 nodes
  
  
  mesh <- inla.mesh.2d(
    mesh0$loc[, 1:2], 
    max.edge=floor(0.388*mean(range.xy)/4), 
    offset=floor(0.259*mean(range.xy)/2),   
    cutoff=floor(floor(0.259*mean(range.xy)/2)/11.5))    
  mesh$n  
  if(plot == TRUE){
    plot(mesh, asp=1, lwd=2)
    plot(spdf, add=TRUE, border=gray(.7))
  } else {
    
  }
  
  return(mesh)
  
}



