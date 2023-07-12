## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
##
##  ClimateData 3.0
##  R Pipelines for climate data request and processing
##
## ------------------------------------
##
##  theMarineDataScientist
##  github.com/jorgeassis
##  medium.com/themarinedatascientist
##  medium.com/@jorgemfa
##
## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
 
nCores <- detectCores() / 2
idwPower <- 2
interpNMax <- 16

mainFolder <- "." # "/home/jorgeassis/Projects/Bio-ORACLE" # "/media/Bathyscaphe/Bio-ORACLE/" 
tempFolder <- "/home/jorgeassis/Projects/Bio-ORACLE/tempFolder/" # "/media/Bathyscaphe/Bio-ORACLE/tempFolder/"
tempFolder <- paste0(tempFolder,"Rnd_",paste0(sample(LETTERS, 3, TRUE),collapse=""))

rawDataFolder <- c("/media/Bathyscaphe/Raw Data [Environment]/")

projectName <- "BioOracle"
exportFolders <- c("/media/Jellyfish/Bio-ORACLE/Decade/","/media/Bathyscaphe/Bio-ORACLE/")
exportFolder <-exportFolders[1]

mainResolution <- 0.05 # 0.05 // 0.25
bathymetryFolder <- "Dependencies/SpatialData/Rasters"
bathymetryFiles <- c("BathymetryDepthMinRes005.tif","BathymetryDepthMeanRes005.tif","BathymetryDepthMaxRes005.tif") # 025 005
outcomeExtent <- c(-180,180,-90,90)

cat("\n")
cat("\n")
cat("# Raster resolution:",res(raster(paste0(bathymetryFolder,"/",bathymetryFiles[1])))[1],"\n")
cat("# Raster extent:",outcomeExtent)

listVariables <- read.csv("Dependencies/listVariables.csv",sep=";", stringsAsFactors = F)
listModels <- read.csv("Dependencies/listModels.csv",sep=";", stringsAsFactors = F)
listModelsFinal <- read.csv("Dependencies/listModelsFinal.csv",sep=";", stringsAsFactors = F)
listExperiments <- read.csv("Dependencies/listExperiments.csv",sep=";", stringsAsFactors = F)
listDataSources <- read.csv("Dependencies/listDataSources.csv",sep=";", stringsAsFactors = F)
listRealms <- c("ocean","ocnBgchem","ocnBgChem","atmos","seaIce") 

## ------------------------------
## ------------------------------

cat("\n")
cat("# Cores attributed:",nCores,"\n")
cat("\n")
cat("# Reading configuration file // Project",projectName,"\n")

## ------------------------------
## ------------------------------
