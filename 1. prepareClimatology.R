## ---------------------------------------------------------------------
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
## ---------------------------------------------------------------------

cat("# ---------------\n")
cat("# Prepare climatology //",depthInterpolation,"\n")
time.i <- Sys.time()

## ------------------------
## ------------------------

options(warn=-1)

## ------------------

climatologyName <- paste0(variable,"_",experiment,"_",firstToupper(outcomePeriodType),ifelse(is.null(outcomePeriodMName),"",paste0(outcomePeriodMName,"_")),"_",paste0(min(outcomePeriodY.available),"_",ifelse(max(outcomePeriodY.available) != outcomePeriodY[length(outcomePeriodY)] & outcomePeriodType == "decade" & outcomePeriodY[length(outcomePeriodY)] - max(outcomePeriodY.available) <= 10 , outcomePeriodY[length(outcomePeriodY)] , max(outcomePeriodY.available) )))
processTempFolder.i <- paste0(tempFolder,"/mainClimatologies/",climatologyName,ifelse(!is.na(model),paste0("/",model,"/"),"/"))

if( processTempFolder.i != processTempFolder ) {
  
  file.rename(processTempFolder,processTempFolder.i)
  processTempFolder <- processTempFolder.i
  
}

## ------------------

shape <- raster(paste0(bathymetryFolder,"/",bathymetryFiles[1]))
shape <- crop(shape,extent(outcomeExtent))
shape[!is.na(shape)] <- 1

## ------------------

bathymetry.min.m <- raster(paste0(bathymetryFolder,"/",bathymetryFiles[1]))
bathymetry.mean.m <- raster(paste0(bathymetryFolder,"/",bathymetryFiles[2]))
bathymetry.max.m <- raster(paste0(bathymetryFolder,"/",bathymetryFiles[3]))

bathymetry.min.m <- crop(bathymetry.min.m,shape)
bathymetry.mean.m <- crop(bathymetry.mean.m,shape)
bathymetry.max.m <- crop(bathymetry.max.m,shape)

# -------------------------

interpCellsNumber <- Which(!is.na(shape),cells=TRUE)
interpCellsRow <- rowFromCell(shape,interpCellsNumber)
interpCellsColumn <- colFromCell(shape,interpCellsNumber)

interpCells <- data.table( Lon = xyFromCell(shape,interpCellsNumber)[,1] ,
                           Lat = xyFromCell(shape,interpCellsNumber)[,2] ,
                           Row = interpCellsRow ,
                           Col = interpCellsColumn ,
                           Cell = interpCellsNumber ,
                           depthSurf = 0,
                           depthMin = bathymetry.min.m[interpCellsNumber] * (-1),
                           depthMean = bathymetry.mean.m[interpCellsNumber] * (-1),
                           depthMax = bathymetry.max.m[interpCellsNumber] * (-1) )

setkey(interpCells, Lon, Lat)

# -------------------------
# -------------------------

for( pred in outcomePeriodPredictors) {   
  
  fileName <- paste0("climatology",pred,depthInterpolation)
  climatologyBin <- paste0(fileName,".bin")
  climatologyDesc <- paste0(fileName,".desc")
  file.remove( list.files(processTempFolder, full.names = TRUE, pattern = fileName) ) 
  climatologyBm <- big.matrix( nrow = nrow(interpCells)  ,ncol = as.numeric(dim.t) , backingpath=processTempFolder , backingfile = climatologyBin, descriptorfile = climatologyDesc)

}

## -----------

interpCellsBM <- as.matrix(interpCells)
colnames(interpCellsBM) <- NULL
file.remove( list.files(processTempFolder, full.names = TRUE, pattern = "cellsStructure") ) 

climatologyCellsBin <- paste0("cellsStructure",".bin")
climatologyCellsDesc <- paste0("cellsStructure",".desc")
climatologyCellsBm <- as.big.matrix(interpCellsBM , backingpath=processTempFolder , backingfile = climatologyCellsBin, descriptorfile = climatologyCellsDesc)
climatologyCellsDesc <- dget( paste0(processTempFolder,"/",climatologyDesc))

## ------------------ 

rm(interpCellsRow); rm(interpCellsColumn); rm(interpCellsNumber)
closeAllConnections()
gc(reset=TRUE)

cat("# End process // ET:", round(difftime(Sys.time(), time.i, units='hours'), digits = 3) ,"hours\n")

## --------------------------------------------------------
## --------------------------------------------------------