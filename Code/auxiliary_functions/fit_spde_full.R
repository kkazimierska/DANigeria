# Input
# stk.full - full stack for the estimation 
# res.spde - inla result of the spde model fit 
# spde - spde model, the same as for the fit

# Output 
# prediction result of the spde fit 

fit_spde_full <- function(stk.full, res.spde, spde){
  formula.model <- count ~ Intercept + f(i, model=spde)-1
  
  pred.res.spde <- inla(formula.model,
                     family = "poisson",
                     control.family = list(link = "log"),
                     data=inla.stack.data(stk.full), ## supply the full data
                     # link=1` to compute the fitted values with the same link function
                     # as the family specified in the model
                     control.predictor=list(compute=TRUE, link = 1,
                                            A=inla.stack.A(stk.full)), ## full
                     control.mode=list(theta=res.spde$mode$theta, restart=FALSE),
                     # use posterior mode from the spde fit 
                     # instead of refitting, the mode for hyper-param is fixed
                     num.threads="4:1")## use mode already found
  return(pred.res.spde)
}

