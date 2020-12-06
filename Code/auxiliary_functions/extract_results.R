# Input 
# res - inla output from the prediction fit 
# tag - label of the prediction indices 
# gproj - inla.mesh.projector output containing prediction lattice
# plot - T or F 
# sp - spatial polygon, only if plot is T

# Output 
# fitted relative risk mean and standard error data frame

extract_results <- function(res, tag,  grpoj, plot = T, sp = NULL){
  pred.ind <- inla.stack.index(stk.full, tag=tag)$data
  
  rr = data.frame(
    mean = p.res.spde$summary.fitted.val[pred.ind, "mean"],
    sd = p.res.spde$summary.fitted.val[pred.ind,"sd"])
  
  if(plot ==T){
    # identify which points of the new grid are inside NG
    ov <- over(SpatialPoints(gproj$lattice$loc, proj4string= crs(sp)), sp)
    # estimates outside of NG make NA, so they are not visible
    rr$mean[is.na(ov)] <- NA 
    rr$sd[is.na(ov)] <- NA 
    
    # mean
    image.plot(list(x=gproj$x, y=gproj$y,
                    z=matrix(rr$mean, nrow = nxy[1])), asp=1,
               legend.mar=3, main="RR mean")
    plot(mapdf.km, add=TRUE, border=gray(.3))
    points(coo, pch=4, cex=0.5)
    # sd 
    image.plot(list(x=gproj$x, y=gproj$y,
                    z=matrix(rr$sd, nrow = nxy[1])), asp=1,
               legend.mar=3, main="RR sd")
    plot(mapdf.km, add=TRUE, border=gray(.3))
    points(coo, pch=4, cex=0.5)
  }
  return(rr)

}
