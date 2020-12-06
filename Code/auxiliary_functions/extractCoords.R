# Input 
# spdf spatial polygon data frame 
# output
# the coordinates of the borders of the sp 

extractCoords <- function(spdf)
{
  results <- list()
  for(i in 1:length(spdf@polygons[[1]]@Polygons))
  {
    results[[i]] <- spdf@polygons[[1]]@Polygons[[i]]@coords
  }
  results <- Reduce(rbind, results)
  results
}