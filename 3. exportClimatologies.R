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
cat("# Export climatology //",depthInterpolation,"\n")
time.i <- Sys.time()

## ------------------------------
## ------------------------------

if( length(list.files(exportFolderClimatologies, pattern=".desc")) == 0) { stop("Climatologies not found") }

# -------------

shape <- raster(paste0(bathymetryFolder,"/",bathymetryFiles[1]))
shape <- crop(shape,extent(outcomeExtent))
shape[!is.na(shape)] <- 1

## ------------------------------
## Generate .nc File

ncFileName <- paste0("climatology",firstToupper(outcomePeriodType),firstToupper(depthInterpolation),".nc")

if( outcomePeriodType == "decade") {
  period <- outcomePeriodY.available[which(substring(outcomePeriodY.available,4,4) == 0)] 
  period <- as.numeric(as.Date(paste0(period,"-01-01"), "%Y-%m-%d") - as.Date("01-01-1970", "%d-%m-%Y"))
  if( substring(outcomePeriodY.available[length(outcomePeriodY.available)],4,4)=="0" ) { period <- period[-length(period)]}
}

if( outcomePeriodType == "year") {
  period <- as.numeric(as.Date(paste0(outcomePeriodY.available,"-01-01"), "%Y-%m-%d") - as.Date("01-01-1970", "%d-%m-%Y"))
}

if(length(period) != dim.t) { stop("Error :: 381")}

if( file.exists(paste0(exportFolderClimatologies,"/",ncFileName))) { file.remove(paste0(exportFolderClimatologies,"/",ncFileName)) }

nc <- ncCreate(path=paste0(exportFolderClimatologies,"/",ncFileName),
               xmax = extent(shape)[2], xmin = extent(shape)[1], xlength = ncol(shape),
               ymax = extent(shape)[4], ymin = extent(shape)[3], ylength = nrow(shape),
               realm = depthInterpolation,
               experiment = experiment,
               period = period)

close.nc( nc )

## ------------------------------

descFiles <- list.files(exportFolderClimatologies,full.names = TRUE, pattern = ".desc")
for( descFile in descFiles) {
  tx <- readLines(descFile)
  writeLines( gsub("//","/",gsub(processTempFolder, exportFolderClimatologies, tx)) , con=descFile)
}

descFiles <- list.files(exportFolderClimatologies,full.names = TRUE, pattern = ".desc")
for( descFile in descFiles) {
  tx <- readLines(descFile)
  tx <- gsub("\\[|\\]","",tx)
  if( TRUE %in% grepl("temp",tx)) { stop("Error :: 301")}
  if(!TRUE %in% grepl(gsub("//","/",gsub("\\[|\\]","",exportFolderClimatologies)), tx)) { stop("Error :: 302")}
  if(!TRUE %in% grepl(paste0("_",min(outcomePeriodY.available),"_",ifelse(max(outcomePeriodY.available) != outcomePeriodY[length(outcomePeriodY)] & outcomePeriodType == "decade" & outcomePeriodY[length(outcomePeriodY)] - max(outcomePeriodY.available) <= 10 , outcomePeriodY[length(outcomePeriodY)] , max(outcomePeriodY.available) )),tx)) { stop("Error :: 303")}
}

## ------------------------------
## Populate .nc File

nc <- open.nc( paste0(exportFolderClimatologies,"/",ncFileName), write=TRUE )

climatologyCells <- dget( paste0(exportFolderClimatologies,"/","cellsStructure.desc"))
climatologyCells <- attach.big.matrix(climatologyCells)

# Correction of layer

correctionLayer <- 0
  
if( rawDataVarUnits == "K" ) { 
  
  rawDataVarUnits <- "degrees_C" 
  correctionLayer <- -273.15

}

# Load climatologies and fillUp .nc files

for( pred in outcomePeriodPredictors) {   

  assign( paste0("climatology",pred) , dget( paste0(exportFolderClimatologies,"/","climatology",pred,depthInterpolation,".desc")) )
  assign( paste0("climatology",pred) , attach.big.matrix(get(paste0("climatology",pred))) )
  
  dim.t <- ncol(get(paste0("climatology",pred)))
  
  # Populate .nc files # Rotate structure
  
  climatologyArray <- array(dim = c(ncol(shape), nrow(shape), dim.t ))

  for( t in 1:dim.t) {
        
      temArray <- shape
      temArray[climatologyCells[,5]] <- get(paste0("climatology",pred))[,t]
      temArray <- temArray + correctionLayer
      climatologyArray[,,t] <- rotateMatrix(raster::as.matrix(temArray)) # as.matrix(temArray)
      
  }
  
  # Correct to maximum allowed value
  
  if(  pred != "Range" & experiment == "Baseline" ) {
    
    climatologyArray[climatologyArray < listVariables[listVariables$name == variable,"minPossibleValue"] ] <- listVariables[listVariables$name == variable,"minPossibleValue"]
    climatologyArray[climatologyArray > listVariables[listVariables$name == variable,"maxPossibleValue"] ] <- listVariables[listVariables$name == variable,"maxPossibleValue"]
    
  }

  if( pred == "Max" ) { predictorLongName <- paste0("Maximum ", variableLongName); predictorStandardName <- "maximum"; varRangeMax <- range(climatologyArray,na.rm=T) }
  if( pred == "ltMax" ) { predictorLongName <- paste0("Long-term maximum ", variableLongName); predictorStandardName <- "long_term_maximum"; varRangeLtMax <- range(climatologyArray,na.rm=T) }
  if( pred == "Mean" ) { predictorLongName <- paste0("Average ", variableLongName); predictorStandardName <- "mean"; varRangeMean <- range(climatologyArray,na.rm=T) }
  if( pred == "ltMin" ) { predictorLongName <- paste0("Long-term minimum ", variableLongName); predictorStandardName <- "long_term_minimum"; varRangeLtMin <- range(climatologyArray,na.rm=T) }
  if( pred == "Min" ) { predictorLongName <- paste0("Minimum ", variableLongName); predictorStandardName <- "minimum"; varRangeMin <- range(climatologyArray,na.rm=T) }
  if( pred == "Range" ) { predictorLongName <- paste0("Range ", variableLongName); predictorStandardName <- "range" }
  
  ncAddVar(nc, varname = paste0(variable,"_",tolower(pred)), varUnits = metadataUnits, long_name = predictorLongName, "crs", compression=9)
  var.put.nc(nc, variable = paste0(variable,"_",tolower(pred)), data = climatologyArray, start = c( 1, 1, 1), count = c( dim(climatologyArray)[1], dim(climatologyArray)[2], dim(climatologyArray)[3]) )

}

close.nc(nc)

## --------------------------

rawMetadata <- read.csv(file=paste0(exportFolderClimatologies,"/rawDataMetadata.txt"),dec=".",sep=";",header=TRUE)

metadata <- data.frame(   descriptor = c("variableLongName" ,
                                         "variableStandardName",
                                         "variable",
                                         "realm",
                                         "units",
                                         "absoluteValues",
                                         "dataSource",
                                         "dataProduct",
                                         "dataDataset",
                                         "experiment",
                                         "ensemble",
                                         "climatologyType",
                                         "climatologyPeriod",
                                         "climatologyPeriodN",
                                         "vectorValues",
                                         "minValues",
                                         "maxValues",
                                         "resolution",
                                         "nCells"),
                          
                          value = c(variableLongName,
                                    variableStandardName,
                                    variable,
                                    ifelse(realm == "Benthic",paste0(realm,firstToupper(depthInterpolation)),realm),
                                    rawMetadata$value[rawMetadata$descriptor == "units"],
                                    listVariables[listVariables$name == variable,4],
                                    rawMetadata$value[rawMetadata$descriptor == "dataSource"],
                                    rawMetadata$value[rawMetadata$descriptor == "dataProduct"],
                                    rawMetadata$value[rawMetadata$descriptor == "dataDataset"],
                                    experiment,
                                    "None",
                                    firstToupper(outcomePeriodType),
                                    paste0(min(outcomePeriodY.available),"-",max(outcomePeriodY.available)),
                                    length(outcomePeriodY.available),
                                    paste0(c("min","ltmin","mean","ltmax","max"), collapse = ", "),
                                    paste0(c( varRangeMin[1], ifelse(exists("varRangeLtMin"),varRangeLtMin[1],NA), varRangeMean[1], ifelse(exists("varRangeLtMax"),varRangeLtMax[1],NA), varRangeMax[1] ), collapse = ", "),
                                    paste0(c( varRangeMin[2], ifelse(exists("varRangeLtMin"),varRangeLtMin[2],NA), varRangeMean[2], ifelse(exists("varRangeLtMax"),varRangeLtMax[2],NA), varRangeMax[2] ), collapse = ", "),
                                    res(shape)[1],
                                    ncell(shape)
                          ) )

write.table(metadata,file=paste0(exportFolderClimatologies,"/",gsub(".nc","Metadata.txt",ncFileName)),col.names = TRUE,sep=";",dec=".",row.names = FALSE,quote=FALSE)

## --------------------------
## --------------------------

cat("# End process // ET:", difftime(Sys.time(), time.i, units='hours') ,"hours\n")

## --------------------------------------------------------
## --------------------------------------------------------