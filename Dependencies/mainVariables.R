## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
##
##                                      #####
##                                ####  #####
##                                ####       
##          ####                         
##         ##################             
##           ##################           
##       #######################
##   ##################################   
##  #######################################
##  ######################################
##  ###################################### 
##  ####################################
##  ##################################     
##  ####################                   
##  ###################                    
##  ##################                     
##  #################                      
##  ###############                                     
##      
##  theMarineDataScientist
##
##  github.com/jorgeassis
##  medium.com/themarinedatascientist
##  medium.com/@jorgemfa
##
## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------
##
##  bioORACLE 3.0
##  R Pipelines for climate data request and processing
##
## -------------------------------------------------------------------------------
## -------------------------------------------------------------------------------

mainGrid <- getGrid(realm,maxDepth)

variableLongName <- listVariables[listVariables$name == variable,"longName"]
experimentLongName <- listExperiments[listExperiments$name == experiment,"longName"]

if(realm == "surface") { interpN <- round(interpNMax / 2) } 
if(realm != "surface") { interpN <- interpNMax } 

# files and folder

if( variable == "mlotst" & experiment == "present") {
  
  rawDataFolder <- paste0(mainRawDataFolder,"/RawDataMyOcean/GLOBAL_REANALYSIS_PHY_001_030")
  rawDataFiles <- list.files(rawDataFolder,full.names = TRUE,recursive=TRUE,pattern = "nc")
  
  rawData <- nc_open( rawDataFiles[1], readunlim=FALSE )
  cat("File:", ncatt_get(rawData,"time")$long_name,"\n")
  cat("Convertion: 1950-01-01")
  cat("\n")
  nc_close( rawData )
  
  filesToKeep <- numeric(0)
  for( f in 1:length(rawDataFiles) ) {
    
    rawData <- nc_open( rawDataFiles[f], readunlim=FALSE )
    rawDataVartime <- ncvar_get( rawData, "time" ) / 24
    rawDataVartime <- as.numeric(format(as.Date(rawDataVartime, origin = "1950-01-01"),"%Y"))
    nc_close( rawData )
    
    if( rawDataVartime %in% outcomePeriodY) { filesToKeep <- c(filesToKeep,f)}
    
  }

  rawDataFiles <- rawDataFiles[filesToKeep]

}

if( variable == "parbottom" & experiment == "present") {

  rawDataFolder <- paste0(mainRawDataFolder,"/RawDataLightAtBottom")
  rawDataFiles <- list.files(rawDataFolder,full.names = TRUE,recursive=TRUE,pattern = "max_00.nc")
  
}


if( variable == "exposureW") {
  
  rawDataFolder <- paste0(mainRawDataFolder,"/RawDataCoastalExposure")
  rawDataFiles <- list.files(rawDataFolder,full.names = TRUE,recursive=TRUE,pattern = "W_map.csv")
  
}

if( variable == "exposureWX") {
  
  rawDataFolder <- paste0(mainRawDataFolder,"/RawDataCoastalExposure")
  rawDataFiles <- list.files(rawDataFolder,full.names = TRUE,recursive=TRUE,pattern = "WXSD_map.csv")
  
}

## ------------------------------
## ------------------------------

