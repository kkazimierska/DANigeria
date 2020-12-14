# Input 
# result - INLA spde fit 
# N - number of samples, default is 2e3 
# Ap - projector matrix for the prediction 
# plot - plot the summary statistics, mean, sd, probability of exceedance 
# sp - spatial polygon data frame, if plot == T
# gproj - inla.mesh.projector output containing the prediction lattice, if plot == T
# coo - coordinates of the centroids of the polygons 

# Output 
# theta.samp - monte carlo samples of the fitted values 
# rr - sumarry statistics, mean, sd and probability of excedance 

MC_samples <- function(result, N = 2e3, Ap, plot = F, sp = NULL, gproj = NULL, coo = NULL){
  
  # control.compute$config = T has to be computed in the result
  samples <- inla.posterior.sample(n=N, result=result, add.names=FALSE)
  # every sample is as a list
  
  # identify intercept and the spatial index 
  # we collect all the names of the latent GMRF
  # first is response, than there is a predictor related to spde, mesh nodes of the random effect u
  # the last is Intercept
  xnames <- rownames(samples[[1]]$latent)
  idx <- lapply(c('Intercept', 'i'), function(nam) ## for each effect
    which(substr(xnames, 1, nchar(nam))==nam)) ## find the index
  
  # indexes are used to collect the the latent field and organize it into a matrix
  # go over the samples and take only the samples of latent field elements: 
  # Intercept and spatial field
  mat.samples <- sapply(samples, function(spl)
    c(Intercept=spl$latent[idx[[1]]],
      u=spl$latent[idx[[2]]]))
  
  # compute the linear predictor usign projection matrix for the prediction 
  eta.g.samp <- as.matrix(cbind(b0=1, s = Ap)%*%mat.samples)
  theta.g.samp <- exp(eta.g.samp)
  
  rr <- list(mean = rowMeans(theta.g.samp),
             sd = matrixStats::rowSds(theta.g.samp),
             p1 = rowMeans(theta.g.samp > 1))
  if(plot ==T){
    ov <- over(SpatialPoints(gproj$lattice$loc, proj4string= crs(sp)), sp)
    # estimates outside of NG make NA, so they are not visible
    for (s in c(1,2,3)) {
      rr[[s]][is.na(ov)] <- NA  # outside the polygon of interest make it NA
      fields::image.plot(list(x=gproj$x, y=gproj$y,
                      z=matrix(rr[[s]], nrow = gproj$lattice$dims[1])), asp=1,
                 legend.mar=3, main=names(rr)[s])
      plot(sp, add=TRUE, border=gray(.3))
      points(coo, pch=4, cex=0.5)
    }
   }
  
  return(list(theta.samp = theta.g.samp, rr = rr))
}
