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

availableFiles <- list.files(path = mainDataFolder , pattern = ".nc" , full.names = TRUE , recursive = TRUE)
availableFilesNames <- list.files(path = mainDataFolder , pattern = ".nc" , full.names = FALSE , recursive = TRUE)
length(availableFiles)

availableFiles <- availableFiles[grepl("/Historical/",availableFiles)]
availableFilesNames <- availableFilesNames[grepl("/Historical/",availableFilesNames)]
length(availableFiles)
length(availableFilesNames)

availableData <- data.frame( variable = rep("",length(availableFiles)) , 
                             model = rep("",length(availableFiles)) , 
                             timeFrame = rep("",length(availableFiles)) , 
                             period = rep("",length(availableFiles)) , 
                             ensemble = rep("",length(availableFiles)) , 
                             experiment = rep("",length(availableFiles)) , 
                             realm3D = rep(FALSE,length(availableFiles)) , 
                             rangeMin = rep(NA,length(availableFiles)) , 
                             rangeMax = rep(NA,length(availableFiles)) , 
                             units = rep("",length(availableFiles)) , 
                             stringsAsFactors = F) 

## ------------------------------
## ------------------------------

for(f in 1:length(availableFiles)) {
  
  divider <- as.numeric(unlist(gregexpr("/",availableFilesNames[f])))
  experiment <- substring(availableFilesNames[f],divider[1]+1,divider[2]-1)
  
  if( ! grepl("RawDataCMIP5",availableFilesNames[f]) ) { next }
  
  dividers <- unlist(gregexpr("_",availableFilesNames[f]))
  variable <- substring(availableFilesNames[f],divider[2]+1,dividers[1]-1)
  timeFrame <- substring(availableFilesNames[f],dividers[1]+1,dividers[2]-1)
  model <- substring(availableFilesNames[f],dividers[2]+1,dividers[3]-1)
  ensemble <- substring(availableFilesNames[f],dividers[4]+1,dividers[5]-1)
  
  divider <- as.numeric(gregexpr(".nc",availableFilesNames[f])[[1]])
  period <- substring(availableFilesNames[f],dividers[5]+1,divider-1)
  
  rawData <- nc_open( availableFiles[f], readunlim=FALSE )
  
  rawDataVarUnits <- ncatt_get(rawData,variable,"units")$value

  get.dim.names(rawData)

  availableData[f,1] <- variable
  availableData[f,2] <- model
  availableData[f,3] <- timeFrame
  availableData[f,4] <- period
  availableData[f,5] <- ensemble
  availableData[f,6] <- experiment
  
  if(is.null(lon.var.name)) { next }  
  
  rawDataLongitude <- ncvar_get( rawData, lon.var.name )
  rawDataLatitude <- ncvar_get( rawData, lat.var.name )
  
  dim.i <- dim(rawDataLongitude)[1]
  dim.j <- dim(rawDataLongitude)[2]
  
  if(is.na(dim.j)) { dim.j <<- dim(rawDataLatitude) }
  
  tryCatch( test.dims <- ncvar_get(rawData,variable, start=c(1,1,1,1), count=c(dim.i,dim.j,1,1)) , error=function(e){ is3DRawData <<- FALSE } )
  tryCatch( test.dims <- ncvar_get(rawData,variable, start=c(1,1,1), count=c(dim.i,dim.j,1)) , error=function(e){ is3DRawData <<- TRUE } )
  
  if( ! is3DRawData ) { dim.z <<- 0 }

  availableData[f,7] <- is3DRawData
  availableData[f,8] <- range(test.dims,na.rm=T)[1]
  availableData[f,9] <- range(test.dims,na.rm=T)[2]
  availableData[f,10] <- rawDataVarUnits
  
  nc_close(rawData)
}

availableData <- availableData[availableData$variable != "",]
View(availableData[sort(availableData$variable,index.return=T)$ix,])
write.csv(availableData,file="Resources/availableData.csv")
