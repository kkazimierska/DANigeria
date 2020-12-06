# Input 
# spde - spde model 
# stack - stack 

# Ouput
# fitted INLA model

fit_spde <- function(spde, stack){
  
  formula.model <- count ~ Intercept + f(i, model=spde)-1
  ##+ f(ic, copy='i', fixed=FALSE)
  
  res.spde <- inla(
    formula.model,
    family = "poisson", 
    data=inla.stack.data(stack),
    E = inla.stack.data(stack)$E, 
    control.predictor=list(
      A=inla.stack.A(stack), compute=TRUE),
    control.compute = list(dic = T, waic = T, cpo = T, config = TRUE),
    control.inla=list(int.strategy='eb'), # eb strategy is to run faster or adaptive 
    num.threads="4:1")
  
}