## ------------------------------------------------------------------ ##
##                                                                    ##
##  ClimateData v3.0                                                  ##
##  R Pipelines for climate data request and processing               ##
##                                                                    ##
## --------------------------------------------                       ##
##                                                                    ##
##  Jorge Assis [ jorgemfa@gmail.com ]                                ##
##  biodiversityDS                                                    ##
##  biodiversityDataScience.com                                       ##
##                                                                    ##
## ------------------------------------------------------------------ ##

cat("# ---------------\n")
cat("# Export rasters //",depthInterpolation,"\n")
time.i <- Sys.time()

## ------------------------------
## ------------------------------

if( exists("availableClimatologies") & experiment != "Baseline" ) { exportFolderClimatologies <- paste0(exportFolder,"/mainClimatologies/",availableClimatologies) }
if( exists("availableClimatologies") & experiment == "Baseline" ) { exportFolderClimatologies <- paste0(exportFolder,"/mainClimatologies/",availableClimatologies) }

ncFileName <- paste0(exportFolderClimatologies,"/climatology",firstToupper(outcomePeriodType),firstToupper(depthInterpolation),".nc")
rawMetadata <- read.csv(file=gsub("\\.nc","Metadata.txt",ncFileName),dec=".",sep=";",header=TRUE)

nc <- nc_open( ncFileName )
climatologyTime <- ncvar_get( nc, "time")
climatologyTime <- substr(as.Date(climatologyTime,origin = "1970-01-01"),1,4)
nc_close(nc)

if(experiment != "Baseline") { outcomePeriodPredictors <- unique(c(outcomePeriodPredictors,"sd")) }

for(layer in outcomePeriodPredictors) {
  
  for(period.n in 1:length(climatologyTime)) {
    
    periodName <- as.numeric(climatologyTime[period.n])
    
    if(! dir.exists( paste0(exportFolder,"/rasterLayers/",experiment,"/",firstToupper(outcomePeriodType),"/") )) { dir.create( paste0(exportFolder,"/rasterLayers/",experiment,"/",firstToupper(outcomePeriodType),"/") ,recursive = TRUE ) }
    
    rasterFile <- paste0(exportFolder,"/rasterLayers/",experiment,"/",firstToupper(outcomePeriodType),"/",variableLongName," ",ifelse(realm == "Benthic",paste0(realm,firstToupper(depthInterpolation)),realm)," ",firstToupper(layer)," ",ifelse(outcomePeriodType == "decade",paste0(periodName,"-",periodName+10),periodName),".tif")
    
    variable.i <- variable
    if( experiment != "Baseline" & variable == "nppv") { variable.i <- "intpp" }
    
    rasterVar <- ncToRaster(ncFile = ncFileName , layer = paste0(variable.i,"_",tolower(layer)), period = period.n, rasterFile = rasterFile)
    writeRaster(rasterVar,file=rasterFile, format="GTiff", overwrite=TRUE)
    
  }
}

metadataFile <- paste0(exportFolder,"/rasterLayers/",experiment,"/",firstToupper(outcomePeriodType),"/",variableLongName," ",ifelse(realm == "Benthic",paste0(realm,firstToupper(depthInterpolation)),realm)," metadata.txt")
write.table(rawMetadata,file=metadataFile,col.names = TRUE,sep=";",dec=".",row.names = FALSE,quote=FALSE)

## --------------------------
## ------------------ 

cat("# End process // ET:", round(difftime(Sys.time(), time.i, units='hours'), digits = 3) ,"hours\n")

## --------------------------------------------------------
## --------------------------------------------------------