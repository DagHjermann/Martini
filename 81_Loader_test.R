
library(loadeR)
library(dplyr)
library(purrr)

# Example from
# https://github.com/SantanderMetGroup/loadeR/wiki/Dataset-definition-and-loading-local-grid-data

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Loading from a single file
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Slow, skip if it is already downloaded
# download.file("http://meteo.unican.es/work/loadeR/data/Iberia_NCEP.tar.gz", 
#               destfile = "C:/data/temp/Iberia_NCEP.tar.gz")

# Extract files from the tar.gz file
untar("C:/data/temp/Iberia_NCEP.tar.gz", exdir = "C:/data/temp/")

# Define dir
dir <- "C:/data/temp/Iberia_NCEP"

# This new directory contains 6 NetCDF files, each corresponding to a different variable:
list.files(dir, pattern = "\\.nc$")

# Open connection
di <- dataInventory("C:/data/temp/Iberia_NCEP/NCEP_pr.nc")

# Check content
str(di)

# Load data
pr <- loadGridData(dataset = "C:/data/temp/Iberia_NCEP/NCEP_pr.nc", var = "pr")

# This function doesn't exist :-)
plotMeanField(pr)

# Load aggregated data
monthlyPr <- loadGridData(dataset = "C:/data/temp/Iberia_NCEP/NCEP_pr.nc", 
                   var = "pr",
                   aggr.m = "sum")

str(monthlyPr)

t <- 1
dat <- list(
  x = monthlyPr$xyCoords$x,
  y = monthlyPr$xyCoords$y,
  z = t(monthlyPr$Data[t,,])
  )
fields::image.plot(dat)
maps::map("world", add = TRUE)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Loading from a dataset (.ncml)
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# Imagine that we want to load several variables, or that the NetCDF files in the directory are of the same variable but for 
# different time periods and we want to load a large period. In those cases, we would need to repeat the procedure above several times pointing to 
# different files and, in addition, a posterior aggregation of the data would be also necessary. This could be avoided creating virtual datasets, 
# which aggregates files and/or variables, by means of a NcML file. We can create the NcML file with function makeAggregatedDataset, that can handle
# automatically all these situations.

# In order to aggregate the NetCDF files to form a unique dataset (NcML), we can just type:

makeAggregatedDataset(source.dir = dir, ncml.file = "C:/data/temp/Iberia_NCEP/Iberia_NCEP.ncml", verbose = TRUE)

system("cat C:/data/temp/Iberia_NCEP/Iberia_NCEP.ncml")

ncep.local <- "C:/data/temp/Iberia_NCEP/Iberia_NCEP.ncml"
di <- dataInventory(ncep.local)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Loading a grid from the aggregated data set
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

T1000 <- loadGridData(dataset = ncep.local, 
                              var = "T@1000", 
                              lonLim = c(-12, 5), 
                              latLim= c(35,45), 
                              season= 6:8, 
                              years = 1981:2000)

str(T1000)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# MartiniControlRun - load dat from the thredds server
#
# Trying to do the same as Ann Kristin in the MartiniControlRun Jypyter file (using Python)
# Go to http://35.233.64.13/hub/login to try MartiniControlRun
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

filepath = 'http://thredds.met.no/thredds/dodsC/metusers/arildb/MARTINI800_prov_v2.ncml' # The OPENDAP URL

# Connect to database
# Slow (ca 1.7 minutes)
t0 <- Sys.time()
filehandle <- dataInventory(filepath )
Sys.time() - t0

# Variable names
names(filehandle)

# Example 
str(filehandle$temp)
str(filehandle$mask_rho)

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
# Spatial domain
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

filehandle$temp$Dimensions$lon$Values %>% diff() %>% mean()
filehandle$temp$Dimensions$lat$Values %>% diff() %>% mean()

expand.grid(range(filehandle$temp$Dimensions$lon$Values), 
            range(filehandle$temp$Dimensions$lat$Values)) %>% 
  plot()

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
# Get a table of the variables
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

get_summ <- function(x, data = filehandle){
  units <- data[[x]]$Units
  data.frame(Name = x, Description = data[[x]]$Description, Units = ifelse(is.null(units), "NA", units), stringsAsFactors = FALSE)
  }

# testing:
# get_summ("u")
# get_summ("mask_rho")
# map_df(c("u", "ubar"), get_summ)
# map_df(names(filehandle)[1:6], get_summ)

map_df(names(filehandle), get_summ)


#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
# Get values
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

x <- loadeR::

# Load data

temp <- loadGridData(dataset = filepath, 
                     var = "temp",
                     lonLim = c(6, 6.02),
                     latLim = c(56.20, 56.22))
# ERROR MESSAGE:
# [2019-08-21 10:48:41] Opening dataset...
# [2019-08-21 10:48:43] The dataset was successfuly opened
# [2019-08-21 10:48:43] Defining geo-location parameters
# Error in .jcall("RJavaTools", "Ljava/lang/Object;", "invokeMethod", cl,  : 
#                   ucar.ma2.InvalidRangeException: last (-1) must be >= first (0)

temp <- loadGridData(dataset = filepath, 
                     var = "temp",
                     lonLim = c(6.02, 6.00),
                     latLim = c(56.20, 56.22))
# same error

temp <- loadGridData(dataset = filepath, 
                     var = "temp",
                     lonLim = c(6.02, 6.00),
                     latLim = c(56.22, 56.20))
# same error

temp <- loadGridData(dataset = filepath, 
                     var = "temp",
                     lonLim = c(6.00, 6.02),
                     latLim = c(56.22, 56.20))
# same error

temp <- loadGridData(dataset = filepath, var = "temp")
# Error in getVerticalLevelPars(grid, level) : 
# Variable with vertical levels: '@level' following the variable name is required
# Possible values: -0.988095238095238, -0.964285714285714, -0.94047619047619, -0.916666666666667, ...

temp <- loadGridData(dataset = filepath, 
                     var = "temp@-0.988095238095238",
                     lonLim = c(6, 6.02),
                     latLim = c(56.20, 56.22))
# same error as the first one (ucar.ma2.InvalidRangeException)


