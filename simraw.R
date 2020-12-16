## Load the libraries
library(ggplot2)
library(INLA)
library(spdep) # poly2nb 
library(sp)  # over

library("RColorBrewer") # pallete
library(matrixStats) # rowMeans, rowQuantiles, etc
library(data.table) # between CI 
library(knitr) # table 

library(rgeoboundaries) # gb_adm
library(sf)
library(raster) # to read the population object

library(dplyr) # left_join
library(maptools) # unionSpatialPolygon
library(fields) # image 

#library(miceadds) # source.all
library(RandomFields) # GF

# schell script inputs, 16 scenarios, 3 models (iid, bym2, spde)
input = commandArgs(trailingOnly=TRUE)
# scenario.idx = as.integer(input[1])
# model = input[2]
scenario.idx = 13
model = "bym2"
print(scenario.idx)
print(model)

## Folder structure 
## should contain folders Results, Code/auxiliary_functions and Data with pop.RData

## Load data 
load("Data/pop.RData")
# source.all("Code/auxiliary_functions/") # miceadds library
source("Code/auxiliary_functions/areal_mesh_param.R")
source("Code/auxiliary_functions/create_matrix.R")
source("Code/auxiliary_functions/create_mesh.R")
source("Code/auxiliary_functions/create_stack.R")
source("Code/auxiliary_functions/fit_spde.R")
source("Code/auxiliary_functions/extractCoords.R")
source("Code/auxiliary_functions/MC_samples.R")
source("Code/auxiliary_functions/scenario_param.R")


## Scenarios 
map <- st_geometry(gb_adm1("nigeria"))  
n = length(map)
map <- as_Spatial(map, IDs =as.character(1:n)) 
# Transform the CRS to km 
map.km <- spTransform(
  map, CRS('+proj=utm +zone=23 +south +ellps=GRS80 +units=km'))

scenarios = scenario_param(sp = map.km)

weights = "pop.int" 
n.rep = 1 
range = scenarios[[scenario.idx]]$range
sigma2 = scenarios[[scenario.idx]]$sigma2u

beta0 = scenarios[[scenario.idx]]$beta0
adm.lev = scenarios[[scenario.idx]]$adm.lev

## Administrative boundaries 

if(adm.lev == 1){
  map <- st_geometry(gb_adm1("nigeria"))  
} else {
  map <- st_geometry(gb_adm2("nigeria"))  
}

n = length(map)
map <- as_Spatial(map, IDs =as.character(1:n)) 
# Transform the CRS to km 
map.km <- spTransform(
  map, CRS('+proj=utm +zone=23 +south +ellps=GRS80 +units=km'))

## Data

pop.dat.ng <- crop(pop.dat.km, map.km)
# Create a SpatialPolygonsDataFrame 
df = data.frame(lat = coordinates(map.km)[,1], 
                lon = coordinates(map.km)[,2],  
                ID =  seq(1:length(map.km)),
                row.names = 1:length(map.km))
# Polygons ids are ID i
mapdf.km = SpatialPolygonsDataFrame(map.km, df) 

## Simulation of the CG model 

# 1-2. GF at the population raster locations
coo = xyFromCell(pop.dat.km, 1:ncell(pop.dat.km)) # coordinates from all cells of raster r
colnames(coo) = c("coords.x1", "coords.x2")
m.model <- RMmatern(nu = 1, var = sigma2, scale = range/2)
u <-  RFsimulate(m.model, x= coo[,1], y= coo[,2])
r.GF = pop.dat.km # copy the raster but change the values 
values(r.GF) <- u@data$variable1 #  assign the data values to the raster object
# 3. linear predictor
eta = beta0 + values(r.GF) 
# 4. relative risk surface
p.true = exp(eta)/(1+exp(eta)) 
pop = values(pop.dat.km)
p.av.true = sum(pop*p.true)/sum(pop) # average probability
rr.true = p.true/p.av.true # relative risk
r.RRt = pop.dat.km
values(r.RRt) = rr.true
r.RRt <- mask(r.RRt, map.km)
# 5. disease counts based on true p and n
ysim = rpois(dim(coo)[1], lambda = pop*p.true) 
r.y <- pop.dat.km
values(r.y) <- ysim
# 6. pop and counts aggregation
mapdf.km@data$pop <- extract(pop.dat.km, map.km, fun = sum, na.rm = TRUE) 
mapdf.km@data$count <- extract(r.y, map.km, fun = sum, na.rm = TRUE) 
# 7. expected counts
E = mapdf.km@data$pop*sum(mapdf.km@data$count)/sum(mapdf.km@data$pop) 
mapdf.km@data$E <- E 


## Define precision matrix once for all simulations

if(model == "bym2"){
  # create rook (right) adjacency matrix, the simplest, condiering neighbours
  adj.r <- poly2nb(map.km, queen = FALSE) # not a queen
  # convert nb object to graph object compatible with INLA, save it and read it
  nb2INLA(file.path("Data/graph", eval(scenario.idx), eval(model), fsep = ""), adj.r)
  ng.graph = inla.read.graph(file.path("Data/graph", eval(scenario.idx), eval(model), fsep = ""))
  
}

## Define mesh, spde model, projector and prediction matrix once for all simulations

if(model == "spde"){
  # mesh
  param = areal_mesh_param(spdf = map.km)
  mesh = create_mesh(spdf = map.km, max.edge = param$max.edge, cutoff = param$cutoff,
                     offset = param$offset, plot = F)
  # spde model
  spde <- inla.spde2.pcmatern(
    mesh, alpha=2,
    prior.sigma=c(1, 0.01),  # P(sigma > 1) = 0.01, 
    prior.range=c(range, 0.5)) # P(r < range) = 0.5
  
  # projector matrix
  if(weights == "equal" ){
    Aa <- create_matrix(sp = map.km, mesh = mesh, weights = "equal")
  } else {
    Aa <- create_matrix(sp = map.km, mesh = mesh, weights = "pop.int", pop.raster = pop.dat.ng)
  }
  
  # prediction lattice and matrix 
  # coo are the coordinates of the population raster
  Ap = inla.spde.make.A(mesh = mesh, loc = coo)
  
}

## Simualtion study 

# globally defined 
# m.model, coo, r.GF, pop, r.RRt, r.y, mapdf.km@data$pop pop.dat.km, map.km, r.RRm

sim.coverage = function(id) {
  
  ## simulate the data 
  u <-  RFsimulate(m.model, x= coo[,1], y= coo[,2]) #  update the GF values
  values(r.GF) <- u@data$variable1
  # 3. linear predictor
  eta = beta0 + values(r.GF) 
  # 4. relative risk surface
  p.true = exp(eta)/(1+exp(eta)) 
  p.av.true = sum(pop*p.true)/sum(pop) # average probability
  rr.true = p.true/p.av.true # relative risk 
  values(r.RRt) = rr.true # update the simulated relative risk values
  r.RRt <- mask(r.RRt, map.km) # contains NA
  # 5. disease counts based on true p and n
  ysim = rpois(dim(coo)[1], lambda = pop*p.true) 
  values(r.y) <- ysim
  # 6. pop and counts aggregation
  mapdf.km@data$count <- extract(r.y, map.km, fun = sum, na.rm = TRUE) 
  # 7. expected counts
  E = mapdf.km@data$pop*sum(mapdf.km@data$count)/sum(mapdf.km@data$pop) 
  mapdf.km@data$E <- E 
  
  ## INFERENCE and prediction
  if(model == "spde"){
    stack.a <- create_stack(tag = "area", matrix = Aa, spdf = mapdf.km, spde = spde)
    res <- fit_spde(spde = spde, stack = stack.a)
    # prediction and CI 
    rrMC = MC_samples(result = res, N = 2e3, Ap = Ap)
    r.RRm <- r.RRt
    values(r.RRm) <- rrMC$rr$mean
    r.RRm = mask(r.RRm, map.km) # contains NA
    
    r.RRql = r.RRqu = r.RRm
    values(r.RRql) = rrMC$rr$q1
    r.RRql = mask(r.RRql, map.km)
    values(r.RRqu) = rrMC$rr$q9
    r.RRqu = mask(r.RRqu, map.km)
    CI <- cbind(values(r.RRql), values(r.RRqu))
  } else {
    if(model == "bym2"){
      formula.model <- count ~ 1 + f(ID, model = "bym2", graph = ng.graph,
                                     constr = T, scale.model = T,
                                     hyper=list(prec=list(prior="pc.prec", param=c(1, 0.01)), # P(sd >1) = 0.01
                                                phi =list(prior="pc", param=c(0.5, 0.5)))) # P(phi <0.5) = 0.5) 
    } else {
      formula.model <- count ~ 1 + f(ID, model = "iid")
    }
    res = inla(formula = formula.model, 
               family = "poisson",
               data = mapdf.km@data, 
               E =  E, 
               control.compute = list(dic = T, waic = T, cpo = T, config = TRUE),
               control.predictor = list(compute = TRUE),
               num.threads="4:1",
               quantiles = c(0.1, 0.9))
    # prediction and CI 
    mapdf.km@data$RR.m <- res$summary.fitted.values[, "mean"]
    mapdf.km@data$RR.qu <- res$summary.fitted.values[, "0.9quant"]
    mapdf.km@data$RR.ql <- res$summary.fitted.values[, "0.1quant"]
    r.RRm <- rasterize(mapdf.km, pop.dat.km, field = "RR.m", fun = "sum") # has NA
    r.RRql <- rasterize(mapdf.km, pop.dat.km, field = "RR.ql", fun = "sum") # CI
    r.RRqu <- rasterize(mapdf.km, pop.dat.km, field = "RR.qu", fun = "sum") # CI 
    CI <- cbind(values(r.RRql), values(r.RRqu))
  }
  
  ## model performance measures
  # absolute error
  ae = abs(values(r.RRt) - values(r.RRm)) # true - posterior mean 
  # credible interval lenght
  cov.len = CI[,2] - CI[,1]
  # binary coverage vector 
  cov.bin <- c()
  for(idx in 1:length(r.RRt)){
    cov.bin[idx] = between(values(r.RRt)[idx], CI[idx, 1], CI[idx, 2])
  }
  
  # save the perf measures for cells inside spatial polygons
  return(list(cov.bin = cov.bin[!is.na(cov.bin)], 
              cov.len = cov.len[!is.na(cov.len)], 
              ae = ae[!is.na(ae)]))
  
}


## Run the model for n.rep number of datasets
res.sim.raw = mclapply(as.list(1:n.rep), mc.cores=3, FUN = sim.coverage)
m = length(res.sim.raw[[1]]$cov.bin)

## Extract model performance matrices of dimension n.rep x m
cov.bin <- matrix(NA, nrow = n.rep, ncol = m) 
cov.len <- matrix(NA, nrow = n.rep, ncol = m) 
ae <- matrix(NA, nrow = n.rep, ncol = m) 

for(id in 1:n.rep){
  # index of the list corresponds to the sim index
  cov.bin[id, ] = res.sim.raw[[id]]$cov.bin
  cov.len[id, ] = res.sim.raw[[id]]$cov.len
  ae[id, ] = res.sim.raw[[id]]$ae
}



## Save the results 
save(cov.bin, cov.len, ae, file = file.path( "Results/mp", eval(scenario.idx), 
                                             eval(model), ".RData", fsep = "")) 
