# Input 
# sp - a spatial polygon data frame, to overlay sp with the mesh nodes, to match polygon id to mesh node 
# mesh - a mesh to create a matrix based on mesh nodes and to which polygons those mesh nodes belongs 
# weights - "equal" or "pop.int", the weights of the matrix for the area are equal or depend on population intensity 
# pop.raster - population raster if weights = "pop.int"


create_matrix <- function(sp, mesh, weights, pop.raster = NULL){
  
  # create Spatial Point object 
  mesh.spt <- SpatialPoints(mesh$loc[, 1:2], sp@proj4string)
  # for each mesh node assign spatial polygon
  loc <- over(x = mesh.spt, y = sp)
  # mesh indexes from 1:(number of mesh nodes)
  id.in <- as.vector(which(!is.na(loc))) # which are not NA
  locin <- mesh$loc[id.in,1:2]
  
  block <- rep(0, nrow(locin))
  for(i in 1:length(sp)){
    block[as.vector(which(!is.na(over(mesh.spt[id.in], sp[i]))))] <- i
  }

  # for all mesh nodes that are not NAs, take the coordinates and which polygon is in
  mesh.df <- data.frame(mlong=locin[,1],
                        mlat=locin[,2],
                        block.id = block)
  
  if(weights == "equal" ){
    # prediction matrix for areal observations 
    Aa <- inla.spde.make.A(
      mesh=mesh, loc=as.matrix(mesh.df[,1:2]),
      block=mesh.df$block.id, block.rescale="sum")
    return(Aa)
  } else {
    # extract pop values at mesh nodes 
    # for the method 'simple' - the values for the cell a point falls in are returned.
    pop.at.meshX <- raster::extract(pop.raster, mesh.spt@coords)
    pop.at.meshX[is.na(pop.at.meshX)] <- 0
    
    # create a data frame of the mesh coordinates and it's population values  
    pop.at.mesh <- data.frame(mesh.spt@coords)
    pop.at.mesh$pop <- pop.at.meshX
    names(pop.at.mesh)[1:2] <- c("mlong", "mlat")
    pop.at.mesh <- plyr::join(mesh.df, pop.at.mesh[id.in,], type="right")
    
    # prediction matrix for the areal observations
    D <- inla.spde.make.A(mesh=mesh, loc=as.matrix(pop.at.mesh[,1:2]), 
                          block=pop.at.mesh$block.id, weights=pop.at.mesh$pop)
    # scaling of the matrix 
    D.tmp <- list()
    D.tmp$D <- D/rowSums(D)
    D.tmp$mesh.weights <- colSums(D)
    D <- D.tmp$D
    return(D)
  }
}





