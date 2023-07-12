## ------------------------------------------------------------------ ##
##                                                                    ##
##  ClimateData v3.0                                                  ##
##  R Pipelines for climate data request and processing               ##
##                                                                    ##
## --------------------------------------------                       ##
##                                                                    ##
##  Jorge Assis [ jorgemfa@gmail.com ]                                ##
##  biodiversityDS                                                    ##
##  biodiversitydatascience.com                                       ##
##                                                                    ##
## ------------------------------------------------------------------ ##

cat("# ---------------\n")
cat("# Ensemble climatologies //",depthInterpolation,"\n")
time.i <- Sys.time()

## ------------------------------
## ------------------------------

## Get matching timeSteps and test for range of values

minValuesListClimatology <- numeric(0)
maxValuesListClimatology <- numeric(0)

for( i in 1:length(availableClimatologiesModels) ) {
  
  metadata <- list.files(availableClimatologiesModels[i],full.names=TRUE,pattern=".txt")
  metadata <- metadata[grepl(firstToupper(depthInterpolation),metadata)]
  metadata <- read.csv(metadata, sep=";")
  timeSteps.i <- metadata$value[ metadata$descriptor == "climatologyPeriod"]
  timeSteps.i <- substr(timeSteps.i,1,4):substr(timeSteps.i,6,9)
  
  if( i == 1) { timeSteps <- timeSteps.i }
  if( i != 1) { timeSteps <- intersect(timeSteps , timeSteps.i) }
  
  minValuesListClimatology <- c(minValuesListClimatology, min(as.numeric(strsplit(metadata[which(metadata$descriptor == "minValues"),2],", ")[[1]]), na.rm=T) )
  maxValuesListClimatology <- c(maxValuesListClimatology, max(as.numeric(strsplit(metadata[which(metadata$descriptor == "maxValues"),2],", ")[[1]]), na.rm=T) )
  
}

## ------------------------------ ##

metadata <- list.files(mainFolderBaselineClimatology,full.names=TRUE,pattern=".txt")
metadata <- metadata[grepl(firstToupper(depthInterpolation),metadata)]
metadata <- read.csv(metadata, sep=";")

minValuesListBaseline <- min(as.numeric(strsplit(metadata[which(metadata$descriptor == "minValues"),2],", ")[[1]]),na.rm=T)
maxValuesListBaseline <- max(as.numeric(strsplit(metadata[which(metadata$descriptor == "maxValues"),2],", ")[[1]]),na.rm=T)
        
## ------------------------------ ##

metadataUnitsClimatologies <- unique(sapply(availableClimatologiesModels,function(x) { metadata.i <- read.csv(paste0(x,"/rawDataMetadata.txt"), sep=";", header=TRUE) ; metadata.i[which(metadata.i$descriptor == "units"),2] } ))
if( length(metadataUnitsClimatologies) > 1) { stop("Error :: 183")}

metadataUnitsBaseline <- read.csv(paste0(mainFolderBaselineClimatology,"/rawDataMetadata.txt"), sep=";", header=TRUE)
metadataUnitsBaseline <- metadataUnitsBaseline[which(metadataUnitsBaseline$descriptor == "units"),2]

## ------------------------------ ##
# Correct climatology

correctClimatology <- FALSE
correctBaseline <- FALSE

if( metadataUnitsClimatologies != "degC" & metadataUnitsClimatologies != "degrees_C" & metadataUnitsBaseline == "degrees_C" ) { 
  stop("Review :: 281")
}

if( metadataUnitsClimatologies == "0.001" & metadataUnitsBaseline == "1e-3" ) { 
  metadataUnitsClimatologies <- "1e-3"
}

if( metadataUnitsClimatologies == "%" & metadataUnitsBaseline == "(0 - 1)" ) { 
  metadataUnitsClimatologies <- "1"
  metadataUnitsBaseline <- "1"
  climatologyCorrect <- 100
  climatologyCorrectType <- "div"
  correctClimatology <- TRUE
  minValuesListClimatology <- minValuesListClimatology / climatologyCorrect
  maxValuesListClimatology <- maxValuesListClimatology / climatologyCorrect
}

if( metadataUnitsClimatologies == "K" & metadataUnitsBaseline == "K" ) { 
  metadataUnitsClimatologies <- "degC"
  metadataUnitsBaseline <- "degC"
}

if( metadataUnitsClimatologies == "%" & metadataUnitsBaseline == "1" ) { 
  metadataUnitsClimatologies <- "1"
  climatologyCorrect <- 100
  climatologyCorrectType <- "div"
  correctClimatology <- TRUE
  minValuesListClimatology <- minValuesListClimatology / climatologyCorrect
  maxValuesListClimatology <- maxValuesListClimatology / climatologyCorrect
}

if( metadataUnitsClimatologies == "kg m-3" & metadataUnitsBaseline == "mg m-3" ) { 
  metadataUnitsClimatologies <- "mg m-3"
  climatologyCorrect <- 1000
  climatologyCorrectType <- "prod"
  correctClimatology <- TRUE
  minValuesListClimatology <- minValuesListClimatology * climatologyCorrect
  maxValuesListClimatology <- maxValuesListClimatology * climatologyCorrect
}

if( metadataUnitsClimatologies == "mol m-3" & metadataUnitsBaseline == "mmol m-3" ) { 
  metadataUnitsClimatologies <- "mmol m-3"
  climatologyCorrect <- 1000
  climatologyCorrectType <- "prod"
  correctClimatology <- TRUE
  minValuesListClimatology <- minValuesListClimatology * climatologyCorrect
  maxValuesListClimatology <- maxValuesListClimatology * climatologyCorrect
}

if( metadataUnitsClimatologies == "mmol m-3" &  metadataUnitsBaseline == "mol m-3" ) { 
  metadataUnitsBaseline <- "mmol m-3"
  baselineCorrect <- 1000
  baselineCorrectType <- "prod"
  correctBaseline <- TRUE
  minValuesListBaseline <- minValuesListBaseline * climatologyCorrect
  maxValuesListBaseline <- maxValuesListBaseline * climatologyCorrect
}

if( metadataUnitsClimatologies == "mol m-2 s-1" & metadataUnitsBaseline == "mg m-2 day-1" ) { 
  metadataUnitsClimatologies <- "mg m-2 day-1"
  climatologyCorrect <- 12010.700000000053 * 60 * 60 * 24 # mol m-2 s-1 to milligrams of Carbon per day
  climatologyCorrectType <- "prod"
  correctClimatology <- TRUE
  minValuesListClimatology <- minValuesListClimatology * climatologyCorrect
  maxValuesListClimatology <- maxValuesListClimatology * climatologyCorrect
}

if( metadataUnitsClimatologies != metadataUnitsBaseline ) { stop("Error :: 0033") } 

## -------------

metadataUnits <- metadataUnitsClimatologies

## ------------------------------ ##
## ------------------------------ ##
## Generate final ensemble .nc File

shape <- raster(paste0(bathymetryFolder,"/",bathymetryFiles[1]))
shape <- crop(shape,extent(outcomeExtent))
shape[!is.na(shape)] <- 1

ncFileName <- paste0("climatology",firstToupper(outcomePeriodType),firstToupper(depthInterpolation),".nc")

if( outcomePeriodType == "decade") {
  
  period <- timeSteps[which(substring(timeSteps,4,4) == 0)] 
  period <- as.numeric(as.Date(paste0(period,"-01-01"), "%Y-%m-%d") - as.Date("01-01-1970", "%d-%m-%Y"))
  if( substring(timeSteps[length(timeSteps)],4,4)=="0" ) { period <- period[-length(period)]}
  
}

if( outcomePeriodType == "year") {
  
  period <- as.numeric(as.Date(paste0(timeSteps,"-01-01"), "%Y-%m-%d") - as.Date("01-01-1970", "%d-%m-%Y"))
  
}

exportFolderClimatologies <- paste0(exportFolder,"/mainClimatologies/",availableClimatologies)
if( file.exists(paste0(exportFolderClimatologies,"/",ncFileName))) { file.remove(paste0(exportFolderClimatologies,"/",ncFileName)) }

nc <- ncCreate(path=paste0(exportFolderClimatologies,"/",ncFileName),
               xmax = extent(shape)[2], xmin = extent(shape)[1], xlength = ncol(shape),
               ymax = extent(shape)[4], ymin = extent(shape)[3], ylength = nrow(shape),
               realm = depthInterpolation,
               experiment = experiment,
               period = period)

close.nc( nc )

nc <- open.nc( paste0(exportFolderClimatologies,"/",ncFileName), write=TRUE )
ncBaseline <- nc_open( paste0(mainFolderBaselineClimatology,"/",ncFileName))

## --------

minValuesList <- ltMinValuesList <- meanValuesList <- ltMaxValuesList <- maxValuesList <- numeric(0)

for( pred in outcomePeriodPredictors[outcomePeriodPredictors != "sd"] ) {   
  
  for( t in 1:length(period)) {
  
    climatologyArray <- array(dim = c(ncol(shape), nrow(shape), length(availableClimatologiesModels) ))
    
    for( m in 1:length(availableClimatologiesModels)) {
      
      ncFileName.m <- paste0(availableClimatologiesModels[m],"/climatology",firstToupper(outcomePeriodType),firstToupper(depthInterpolation),".nc")
      nc.m <- nc_open( ncFileName.m, readunlim=FALSE )
      climatologyArray[,,m] <- ncvar_get( nc.m, paste0(variable,"_",tolower(pred)), start=c(1,1,t), count=c(ncol(shape),nrow(shape),1))
      nc_close( nc.m )
      
    }
    
    climatologyArray <- arma3DMean(climatologyArray,3)
    
    if( correctClimatology ) {
      if( climatologyCorrectType == "prod" ) {
        climatologyArray <- climatologyArray * climatologyCorrect
      } 
      if( climatologyCorrectType == "div" ) {
        climatologyArray <- climatologyArray / climatologyCorrect
      }
    }
      
    # Apply to baseline
    
    climatologyArrayBaseline.time <- ncvar_get( ncBaseline, "time")
    climatologyArrayBaseline.time <- substr(as.Date(climatologyArrayBaseline.time,origin = "1970-01-01"),1,4)
    climatologyArrayBaseline.time <- which(as.numeric(climatologyArrayBaseline.time) %in% baselinePeriod)
    
    if( length(climatologyArrayBaseline.time) == 1 ) {
      
      climatologyArrayBaseline <- ncvar_get( ncBaseline, paste0(variable.m,"_",tolower(pred)), start=c(1,1,climatologyArrayBaseline.time), count=c(ncol(shape),nrow(shape),1))
      
    }
    
    if( length(climatologyArrayBaseline.time) > 1 ) {
      
      variable.n <- ifelse(variable.m == "nppv" & "PP_mean" %in% names(ncBaseline$var),"PP",variable.m)
      climatologyArrayBaseline <- ncvar_get( ncBaseline, paste0(variable.n,"_",tolower(pred)), start=c(1,1,min(climatologyArrayBaseline.time)), count=c(ncol(shape),nrow(shape),length(climatologyArrayBaseline.time)))
      climatologyArrayBaseline <- arma3DMean(climatologyArrayBaseline,3)
      
    }
    
    if( correctBaseline ) {
      if( baselineCorrectType == "prod" ) {
        climatologyArrayBaseline <- climatologyArrayBaseline * baselineCorrect
      } else { stop("Review :: 445")}
    }

    climatologyArray <- climatologyArrayBaseline + climatologyArray
    
   # if( t == 1 & pred == "Mean" & length(unique(log10Ceiling(c( mean(climatologyArray,na.rm=T) , mean(climatologyArrayBaseline,na.rm=T))))) > 1 ) { stop("Error :: 333") }
   # if( t == 1 & pred == "Mean" & length(unique(log10Ceiling(c( max(climatologyArray,na.rm=T) , max(climatologyArrayBaseline,na.rm=T))))) > 1 ) {  stop("Error :: 334") }
    
    # Correct to maximum allowed value
    
    if(  pred != "Range" ) {
      
      climatologyArray[climatologyArray < listVariables[listVariables$name == variable,"minPossibleValue"] ] <- listVariables[listVariables$name == variable,"minPossibleValue"]
      climatologyArray[climatologyArray > listVariables[listVariables$name == variable,"maxPossibleValue"] ] <- listVariables[listVariables$name == variable,"maxPossibleValue"]
      
    }
    
    if( t == 1) {
      
      if( pred == "Max" ) { predictorLongName <- paste0("Maximum ", variableLongName); predictorStandardName <- "maximum"; varRangeMax <- c(min(climatologyArray,na.rm=T),max(climatologyArray,na.rm=T)) }
      if( pred == "ltMax" ) { predictorLongName <- paste0("Long-term maximum ", variableLongName); predictorStandardName <- "long_term_maximum"; varRangeLtMax <- c(min(climatologyArray,na.rm=T),max(climatologyArray,na.rm=T)) }
      if( pred == "Mean" ) { predictorLongName <- paste0("Average ", variableLongName); predictorStandardName <- "mean"; varRangeMean <- c(min(climatologyArray,na.rm=T),max(climatologyArray,na.rm=T)) }
      if( pred == "ltMin" ) { predictorLongName <- paste0("Long-term minimum ", variableLongName); predictorStandardName <- "long_term_minimum"; varRangeLtMin <- c(min(climatologyArray,na.rm=T),max(climatologyArray,na.rm=T)) }
      if( pred == "Min" ) { predictorLongName <- paste0("Minimum ", variableLongName); predictorStandardName <- "minimum"; varRangeMin <- c(min(climatologyArray,na.rm=T),max(climatologyArray,na.rm=T)) }
      if( pred == "Range" ) { predictorLongName <- paste0("Range ", variableLongName); predictorStandardName <- "range" }
      
      ncAddVar(nc, varname = paste0(variable,"_",tolower(pred)), varUnits = metadataUnitsClimatologies, long_name = predictorLongName, "crs", compression=9)
      
    }
    
    var.put.nc(nc, variable = paste0(variable,"_",tolower(pred)), data = climatologyArray, start = c( 1, 1, t), count=c(ncol(shape),nrow(shape),1))

    if( pred == "Max" ) { maxValuesList <- range(c(maxValuesList,range(climatologyArray,na.rm=T))) }
    if( pred == "ltMax" ) { ltMaxValuesList <- range(c(ltMaxValuesList,range(climatologyArray,na.rm=T))) }
    if( pred == "Mean" ) { meanValuesList <- range(c(meanValuesList,range(climatologyArray,na.rm=T))) }
    if( pred == "ltMin" ) { ltMinValuesList <- range(c(ltMinValuesList,range(climatologyArray,na.rm=T))) }
    if( pred == "Min" ) { minValuesList <- range(c(minValuesList,range(climatologyArray,na.rm=T))) }
    
  }

}

## ------------------------------ ##
# Standard deviation of models

for( t in 1:length(period)) {
  
  climatologyArray <- array(dim = c(ncol(shape), nrow(shape), length(availableClimatologiesModels) ))
  
  for( m in 1:length(availableClimatologiesModels)) {
    
    ncFileName.m <- paste0(availableClimatologiesModels[m],"/climatology",firstToupper(outcomePeriodType),firstToupper(depthInterpolation),".nc")
    nc.m <- nc_open( ncFileName.m, readunlim=FALSE )
    climatologyArray[,,m] <- ncvar_get( nc.m, paste0(variable,"_","mean"), start=c(1,1,t), count=c(ncol(shape),nrow(shape),1))
    nc_close( nc.m )
    
  }
  
  cl <- makeCluster(nCores)
  climatologyArray <- parApply(cl = cl, X = climatologyArray, MARGIN = 1:2, FUN = sd)
  stopCluster(cl)
  closeAllConnections()
  
  if( t == 1) {
    predictorLongName <- paste0("Standard deviation ", variableLongName); predictorStandardName <- "std"; varRangeSD <- c(min(climatologyArray,na.rm=T),max(climatologyArray,na.rm=T))
    ncAddVar(nc, varname = paste0(variable,"_","sd"), varUnits = metadataUnits, long_name = predictorLongName, "crs", compression=9) 
  }
  
  var.put.nc(nc, variable = paste0(variable,"_","sd"), data = climatologyArray, start = c( 1, 1, t), count=c(ncol(shape),nrow(shape),1))
  
}

close.nc(nc)
nc_close( ncBaseline )

## ------------------------------ ##

modelList <- character(0)
for( i in 1:length(availableClimatologiesModels) ) {
  metadata.i <- read.csv(paste0(availableClimatologiesModels[i],"/rawDataMetadata.txt"), sep=";", header=TRUE) 
  modelList <- c(modelList,metadata.i[which(metadata.i$descriptor == "dataSource"),2])
}

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
                                    ifelse( variable == "intpp","nppv",variable),
                                    ifelse(realm == "Benthic",paste0(realm,firstToupper(depthInterpolation)),realm),
                                    metadataUnitsClimatologies,
                                    listVariables[listVariables$name == variable,4],
                                    paste0(modelList,collapse = ", "),
                                    "Ensemble",
                                    "Ensemble",
                                    experiment,
                                    length(modelList),
                                    firstToupper(outcomePeriodType),
                                    paste0(min(as.numeric(substr(as.Date(period, origin = "1970-01-01"),1,4))),"-",max(as.numeric(substr(as.Date(period, origin = "1970-01-01"),1,4)))),
                                    length(period),
                                    paste0(c("min","ltmin","mean","ltmax","max"), collapse = ", "),
                                    paste0(c( minValuesList[1], ltMinValuesList[1], meanValuesList[1], ltMaxValuesList[1], maxValuesList[1] ), collapse = ", "),
                                    paste0(c( minValuesList[2], ltMinValuesList[2], meanValuesList[2], ltMaxValuesList[2], maxValuesList[2] ), collapse = ", "),
                                    res(shape)[1],
                                    ncell(shape)
                          ) )

metadata$value <- gsub(";", ",", metadata$value, fixed = TRUE)
write.table(metadata,file=paste0(exportFolderClimatologies,"/",gsub(".nc","Metadata.txt",ncFileName)),col.names = TRUE,sep=";",dec=".",row.names = FALSE,quote=FALSE)

## ------------------ 

closeAllConnections()
gc(reset=TRUE, full = TRUE)

cat("# End process // ET:", round(difftime(Sys.time(), time.i, units='hours'), digits = 3) ,"hours\n")

## --------------------------------------------------------
## --------------------------------------------------------