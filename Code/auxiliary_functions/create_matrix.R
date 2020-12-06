# Input 
# spdf - a spatial polygon data frame, to overlay sp with the mesh nodes, to match polygon id to mesh node 
# mesh - a mesh to create a matrix based on mesh nodes and to which polygons those mesh nodes belongs 
# weights - "equal" or "pop.int", the weights of the matrix for the area are equal or depend on population intensity 
# pop.raster - population raster if weights = "pop.int"


create_matrix <- function(spdf, mesh, weights, pop.raster = NULL){
  
  # create Spatial Point object 
  mesh.spt <- SpatialPoints(mesh$loc[, 1:2], spdf@proj4string)
 
   # for each node, assign the df from the sp
  loc <- over(x = mesh.spt, y = spdf)
 
   # mesh indexes from 1:(number of mesh nodes)
  id.in <- which(!is.na(loc[,1])) # which are not NA
 
   # for all mesh nodes that are not NAs, assign a polygon ID to which the node falls
   block.id <- pmatch(
    loc$name[id.in], # indices of the geocodes that are not NA
    spdf@data$name,  # spatial polygons
    duplicates.ok=TRUE) # can polygons be used more than once 
  
  # for all mesh nodes that are not NAs, take the coordinates and which polygon is in
  mesh.df <- data.frame(mlong=mesh$loc[id.in, 1],
                        mlat=mesh$loc[id.in, 2],
                        block.id = block.id)
  
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





