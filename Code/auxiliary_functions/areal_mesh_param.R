# Define the max.edge and cutoff for initial mesh based on the boundary
# Define the max.edge and cutoff for the mesh based on the initial mesh locations
# Define the offset for the mesh based on the initial mesh locations

# Input 
# spdf - spatial polygon or spatial polygon data frame

# Output 
# list of vectors max.edge, cutoff and number offset 


areal_mesh_param <- function(spdf){
  
  ext <- extent(coordinates(spdf))
  range.xy = c(ext@xmax-ext@xmin, ext@ymax-ext@ymin)
  max.edge = c(floor(0.078*mean(range.xy)/4), floor(0.388*mean(range.xy)/4)) 
  cutoff = c(floor(0.078*mean(range.xy)/4)/2, floor(floor(0.259*mean(range.xy)/2)/11.5))
  offset = floor(0.259*mean(range.xy)/2)
  
  ls = list(max.edge = max.edge, 
       cutoff = cutoff, 
       offset = offset)
  
  return(ls)
}
