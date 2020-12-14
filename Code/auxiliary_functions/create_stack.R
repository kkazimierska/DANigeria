# Input 
# tag name - label to extract the fitted values 
# matrix - projector matrix for the values to be estimated 
# spdf - spatial polygon data frame with the response and offset values, df names are count and E
# spde - spde model, inla output
# mesh - mesh triangulation, inla output


create_stack <- function(tag, matrix, spdf, spde){
  
  # the A matrix has two elements: A = list(Aa,1)
  # the effect has two elements: first is the set of indexes of and random field multipled by A; 
  # the second is the intercept and this is multiplied by 1, 
  # the matrix Aa has number of column equal to the size of the spatial index and 
  # the number of rows equal to the number of observations
  
  stack.a <- inla.stack(tag=tag, # name obs
                        A=list(matrix,1),
                        data = list(count = spdf@data$count, 
                                    E = spdf@data$E),
                        effects = list(i=1:spde$n.spde,
                                       Intercept = rep(1, nrow(spdf@data))))
  return(stack.a)
}
