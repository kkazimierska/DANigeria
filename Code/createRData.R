library(raster)
lap_data <- "/Users/kazimikk/Dropbox/Paula/Project/Data/"

# disease counts level 1
YF1 <- read.csv(file = file.path(lap_data, "YF1.csv", fsep = ""), sep = ";")

# disease counts level 2
YF2 <- read.csv(file = file.path(lap_data, "YF2.csv", fsep = ""), sep = ";")

pop.dat <- raster(paste0(lap_data, 
                         "wopr/NGA/population/v1.2/NGA_population_v1_2_gridded/NGA_population_v1_2_gridded.tif"))

pop.dat <- setMinMax(pop.dat) # take some time
pop.dat <- aggregate(pop.dat, fun = "sum", fact = 80) 
# 26 k cells
pop.dat.km = projectRaster(pop.dat, 
                           crs = CRS('+proj=utm +zone=23 +south +ellps=GRS80 +units=km'))
pop.dat.km[is.na(values(pop.dat.km))] = 0

save(YF1, YF2, pop.dat.km, file = file.path(lap_data, "ng.RData", fsep = "")) 
