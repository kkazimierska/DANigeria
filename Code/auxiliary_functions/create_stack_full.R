# Input 
# stk.est - stack for the observed data 
# tag - label for the prediction 
# gproj - inla.mesh.projector outpout that contain prediction lattice and projector matrix 
# spde - spde model defined for estimation 
# Output 
# full stack - stack to fit the joint model for estimation and prediction 

create_stack_full <- function(tag, gproj, stk.est, spde){
  # build  a prediction stack similarly as a estimation stack for areal data 
  stack.pred <- inla.stack(tag=tag, # name prediction indices
                           A=list(gproj$proj$A,1),
                           data = list(count = NA, E = NA), ## response as NA
                           effects = list(i=1:spde$n.spde,
                                          Intercept = rep(1, dim(gproj$lattice$loc)[1])))
  
  stk.full <- inla.stack(stk.est, stack.pred) 
  return(stk.full)
}
