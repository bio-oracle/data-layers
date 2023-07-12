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
cat("# Prepare raw data // \n")
time.i <- Sys.time()

## ------------------------------

if( dataSource == "myOcean" ) {
  
  nc.files <- list.files(mainRawDataFolder,full.names = TRUE, recursive=TRUE,pattern = "nc")
  
}

## --------------------

if( dataSource == "cmip6" ) {
  
  nc.files <- list.files(mainRawDataFolder,full.names = TRUE, recursive=TRUE,pattern = "nc")
  nc.files <- nc.files[grepl(model,nc.files)]
  nc.files <- nc.files[grepl(experiment,nc.files)]
  nc.files <- nc.files[grepl(paste0(variable,"_"),nc.files)]

}

## --------------------

if( dataSource == "GlobColour" ) { 
  
  nc.files <- list.files(mainRawDataFolder,full.names = TRUE, recursive=TRUE,pattern = "nc")
  nc.files <- nc.files[grepl(paste0("\\.",variable),nc.files)]
  
}

## --------------------

if( dataSource == "ecmwf" ) { 
  
  nc.files <- list.files(mainRawDataFolder,full.names = TRUE, recursive=TRUE,pattern = "nc")

}

## ------------------------------
## ------------------------------

raw.data <- nc_open( nc.files[1], readunlim=FALSE )
rawDataIs3D <- NULL
singleFile <- FALSE

tryCatch( if(length(ncvar_get( raw.data, variable.i , start=c(1,1,1,1), count=c(1,1,1,1) )) == 1) { rawDataIs3D <<- TRUE } , error=function(e) { errorHand <<- TRUE })
tryCatch( if(length(ncvar_get( raw.data, variable.i , start=c(1,1,1), count=c(1,1,1) )) == 1) { rawDataIs3D <<- FALSE } , error=function(e) { errorHand <<- TRUE })
tryCatch( if(length(ncvar_get( raw.data, variable.i , start=c(1,1), count=c(1,1) )) == 1) { rawDataIs3D <<- FALSE; singleFile <<- TRUE } , error=function(e) { errorHand <<- TRUE })

getDimNames(nc.files[1])
getDimSizes(nc.files[1])

if(  rawDataIs3D ) { raw.data.var <- ncvar_get( raw.data, variable.i , start=c(1,1,1,1), count=c(dim.i,dim.j,dim.z,1) ) }
if( !rawDataIs3D & !singleFile ) { raw.data.var <- ncvar_get( raw.data, variable.i , start=c(1,1,1), count=c(dim.i,dim.j,1) ) }
if( !rawDataIs3D & singleFile ) { raw.data.var <- ncvar_get( raw.data, variable.i , start=c(1,1), count=c(dim.i,dim.j) ) }

if(dimNameDepth == "NA") { depthVarUnits <- dimNameDepth }
if(dimNameDepth != "NA") { depthVarUnits <- getDimUnits(nc.files[1],dimNameDepth) }

rawDataVarUnits <- getDimUnits(nc.files[1],variable.i)
geoVarUnits <- getDimUnits(nc.files[1],dimNameLon)

fillValue <- ncatt_get(raw.data,variable.i,"_FillValue")$value
raw.data.var[raw.data.var == fillValue] <- NA

if( ! grepl("degrees",geoVarUnits) ) { stop("Error :: 995")}

## --------------------

## If ice, fill up the ocean space with zeros

correct.raw.data.matrix <- NULL

if( sum(raw.data.var == 0, na.rm=T) < sum(raw.data.var != 0, na.rm=T) & (variable == "siconc" | variable == "sithick" | variable == "sicov" | variable == "sic" | variable == "sit") ) {
  
  if( dataSource == "myOcean" ) {
    
    correct.raw.data.matrix <- ncvar_get( raw.data, "thetao" , start=c(1,1,1,1), count=c(dim.i,dim.j,1,1) )
    correct.raw.data.matrix[correct.raw.data.matrix != 0] <- 1
    
  }
  
  if( dataSource == "cmip6" ) {
    
    files.to.screen <- list.files(rawDataFolder, recursive=TRUE, pattern=".nc", full.names = TRUE)
    files.to.screen <- files.to.screen[grepl(model,files.to.screen)]
    files.to.screen <- files.to.screen[grepl(experiment,files.to.screen)]
    files.to.screen <- files.to.screen[grepl("thetao",files.to.screen)][1]
    files.to.screen <- nc_open( files.to.screen, readunlim=FALSE )
    
    correct.raw.data.matrix <- ncvar_get( files.to.screen, "thetao", start=c(1,1,1,1), count=c(dim(raw.data.var)[1],dim(raw.data.var)[2],1,1))
    correct.raw.data.matrix[correct.raw.data.matrix != 0] <- 1
    
    if(dim(correct.raw.data.matrix)[1] != dim(raw.data.var)[1] | dim(correct.raw.data.matrix)[2] != dim(raw.data.var)[2]) { stop("Error :: 499 ") }
    
    nc_close( files.to.screen )
    
  }
  
  for( ii in 1:nrow(correct.raw.data.matrix)) {
    for( jj in 1:ncol(correct.raw.data.matrix)) {
      if( ! is.na(correct.raw.data.matrix[ii,jj]) & is.na(raw.data.var[ii,jj]) ) { raw.data.var[ii,jj] <- 0 }
    }
  }
  
}
  
## --------------------

if( rawDataVarUnits == "degC" ) { rawDataVarUnits <- "degrees_C" }

## --------------------

metadataUnits <- rawDataVarUnits
metadataProduct <- ncatt_get( raw.data , 0)$parent_source_id
metadataDataset <- ncatt_get( raw.data , 0)$source
metadataSource <- ncatt_get( raw.data , 0)$source_id

if( is.null(metadataProduct)) { metadataProduct <- dataSource }
if( is.null(metadataDataset)) { metadataDataset <- dataSource }
if( is.null(metadataSource)) { metadataSource <- dataSource }

errorCatch <- FALSE
tryCatch( time <- nc.get.time.series(raw.data), error = function(e){ errorCatch <<- TRUE } )    
if( ! errorCatch ) {  time <- as.Date(as.character(time)) }
if( errorCatch ) { time <- as.Date(ncvar_get( raw.data, "time"),origin="1900-01-01") }

## --------------------

# if( ! is.na(dim.z) & dim.z != length(time) ) { rawDataIs3D <- TRUE } else { rawDataIs3D <- FALSE }

## --------------------

# if( variable == "siconc" | variable == "sithick" | variable == "fice" | variable == "hice" ) { stop("Review :: 021") } 

## --------------------

raw.data.longitude <- ncvar_get( raw.data, dimNameLon )
raw.data.latitude <- ncvar_get( raw.data, dimNameLat )

if( ! rawDataIs3D | (length(depthsToInterpolate) == 1 & "depthSurf" %in% depthsToInterpolate) ){ raw.data.depth <- 0 }

if(   rawDataIs3D & length(depthsToInterpolate) != 1 ){ 
  raw.data.depth <- ncvar_get( raw.data, dimNameDepth ) 
  if( depthVarUnits == "centimeters" ){ raw.data.depth <- raw.data.depth / 100 }
  if(raw.data.depth[1] > 10) { stop("Error: 029")}
  if( depthVarUnits != "centimeters" & depthVarUnits != "m" & depthVarUnits != "0") { stop("Error: 024")}
  
}

nc_close( raw.data )

if( sum(raw.data.longitude[raw.data.longitude < 9999] > 180, na.rm=T) > 0 ) { raw.data.longitude[which(raw.data.longitude > 180)] <- raw.data.longitude[which(raw.data.longitude > 180)] - 360 }
if( min(raw.data.longitude[raw.data.longitude < 9999], na.rm=T) < -299 & max(raw.data.longitude[raw.data.longitude < 9999], na.rm=T) < 60 ) { raw.data.longitude[which(raw.data.longitude < -180)] <- raw.data.longitude[which(raw.data.longitude < -180)] + 360 }
if( min(raw.data.longitude[raw.data.longitude < 9999], na.rm=T) < -180 | max(raw.data.longitude[raw.data.longitude < 9999], na.rm=T) > 180 | min(raw.data.latitude[raw.data.latitude < 9999], na.rm=T) < -90 | max(raw.data.latitude[raw.data.latitude < 9999], na.rm=T) > 90  ) { stop("Error :: 810") }

raw.data.longitude[raw.data.longitude > 9999] <- NA
raw.data.latitude[raw.data.latitude > 9999] <- NA

## --------------------

timeListAvailable <- data.frame()

for( f in 1:length(nc.files)) {
  
  raw.data <- nc_open( nc.files[f], readunlim=FALSE )
  
  errorCatch <- FALSE
  tryCatch( time <- nc.get.time.series(raw.data), error = function(e){ errorCatch <<- TRUE } )    
  if( ! errorCatch ) {  time <- as.Date(as.character(time)) }
  if( errorCatch ) { time <- as.Date(ncvar_get( raw.data, "time"),origin="1900-01-01") }
  nc_close( raw.data )
  
  if(is.na(time) & singleFile ) { 
    
    time <- substr(nc.files[f] , unlist(gregexpr(variable.i,nc.files[f])) + nchar(variable.i) + 1 , unlist(gregexpr("-",nc.files[f])) - 1 ) 
    time <- as.Date(paste0(substr(time,1,4),"-",str_pad(substr(time,5,6), 2, pad = "0"),"-",str_pad(substr(time,8,9), 2, pad = "0")))
    
  }
  
  timeListAvailable <- rbind(timeListAvailable,data.frame(file=nc.files[f],period=time))
  
}

timeListAvailable <- timeListAvailable[sort(as.numeric(timeListAvailable$period), index.return=TRUE)$ix,]
timeListAvailable <- timeListAvailable[which(! duplicated(timeListAvailable[,2])),]

## ----------------

if( experiment == "Baseline" ) {
  timeListAvailable <- timeListAvailable[which(as.numeric(substr(timeListAvailable$period,1,4)) %in% outcomePeriodY), ]
  timeListAvailable <- timeListAvailable[which(as.numeric(substr(timeListAvailable$period,6,7)) %in% outcomePeriodM), ]
}

if( experiment != "Baseline" ) {
  timeListAvailable <- timeListAvailable[which(as.numeric(substr(timeListAvailable$period,1,4)) %in% unique(c(baselinePeriod,outcomePeriodY))), ]
  timeListAvailable <- timeListAvailable[which(as.numeric(substr(timeListAvailable$period,6,7)) %in% outcomePeriodM), ]
}

timeListAvailable <- timeListAvailable[1:which(as.numeric(substr(timeListAvailable$period,6,7)) == 12)[length(which(as.numeric(substr(timeListAvailable$period,6,7)) == 12))],]
timeListAvailableY <- unique(substr(timeListAvailable$period,1,4))
timeListAvailableM <- unique(substr(timeListAvailable$period,6,7))
  
## --------------------

cat("#",paste0(min(as.numeric(timeListAvailableY)),"-",max(as.numeric(timeListAvailableY))),"available for outcome period",paste0(min(as.numeric(outcomePeriodY)),"-",max(as.numeric(outcomePeriodY))),"\n")
outcomePeriodY.n <- which(outcomePeriodY %in% timeListAvailableY)
outcomePeriodY.available <- outcomePeriodY[outcomePeriodY.n]

dim.t <- length(outcomePeriodY.available)
if( dim.t == 0 ) { stop("Error :: 419 ") }

## --------------------
# Decadal time steps

if(outcomePeriodType == "decade" ) {
  
    climatologyPeriod <- outcomePeriodY.available
    
    if( substring(climatologyPeriod[length(climatologyPeriod)],4,4) != 0) { 
      decadalTimeSteps <- seq(climatologyPeriod[1],climatologyPeriod[length(climatologyPeriod)],by=10)
      
      if(decadalTimeSteps[length(decadalTimeSteps)] < climatologyPeriod[length(climatologyPeriod)] ) { decadalTimeSteps <- c(decadalTimeSteps,decadalTimeSteps[length(decadalTimeSteps)] + 10)}
      
      decadalTimeSteps <- data.frame(from = decadalTimeSteps[-length(decadalTimeSteps)], to = decadalTimeSteps[-1] - 1 )
      decadalTimeSteps[nrow(decadalTimeSteps),2] <- climatologyPeriod[length(climatologyPeriod)]
      decadalTimeSteps <- matrix(as.matrix(decadalTimeSteps),ncol=2)
      
    }
    
    if( substring(climatologyPeriod[length(outcomePeriodY.available)],4,4) == 0) { 
      decadalTimeSteps <- seq(climatologyPeriod[1],climatologyPeriod[length(climatologyPeriod)],by=10)
      decadalTimeSteps <- data.frame(from = decadalTimeSteps[-length(decadalTimeSteps)], to = decadalTimeSteps[-1] - 1 )
      decadalTimeSteps <- matrix(as.matrix(decadalTimeSteps),ncol=2)
    }
    if( length(which(apply(decadalTimeSteps,1,diff) == 0)) > 0) {
      decadalTimeSteps <- decadalTimeSteps[-which(apply(decadalTimeSteps,1,diff) == 0),]
    }
    
    timeListAvailableDecDaysSince <- as.numeric(as.Date(paste0(decadalTimeSteps[,1],"-01-01"), "%Y-%m-%d"))
    dim.t <- length(timeListAvailableDecDaysSince)
    
}

## --------------------------------------------------------
## --------------------------------------------------------

if( ! rawDataIs3D | is.na(dim(raw.data.var)[3]) ) { raw.data.var <- array(data=raw.data.var,dim=c(nrow(raw.data.var),ncol(raw.data.var),1))  }
if( dataSource == "GlobColour" | ! useArmadillo ) { interpSurfaceMetric <- data.table(as.data.table(raw.data.var , na.rm=FALSE)) }
if( dataSource != "GlobColour" & useArmadillo ) { interpSurfaceMetric <- data.table(as.data.table(raw.data.var , na.rm=TRUE)) }

## -----------

if( rawDataIs3D ) { 
  
  rawDataVarTest <- interpSurfaceMetric[ which(interpSurfaceMetric[,3] == 1 ),"value"][[1]]
  interpSurfaceMetric <- data.frame(Var1=interpSurfaceMetric[,1],Var2=interpSurfaceMetric[,2],Var3=interpSurfaceMetric[,3],Lon=NA,Lat=NA,Depth=NA)
  
} 

if( ! rawDataIs3D ) {
  
  interpSurfaceMetric <- interpSurfaceMetric[ which(interpSurfaceMetric[,3] == 1),]
  rawDataVarTest <- interpSurfaceMetric[,"value"][[1]]
  interpSurfaceMetric <- data.frame(Var1=interpSurfaceMetric[,1],Var2=interpSurfaceMetric[,2],Var3=1,Lon=NA,Lat=NA,Depth=0)
  
}

colnames(interpSurfaceMetric) <- c("Var1","Var2","Var3","Lon","Lat","Depth")

## --------------------------------------------------------

multdimCoordinates <- ! is.na(dim(raw.data.longitude)[2])

if( multdimCoordinates ) {
  
  cl.2 <- makeCluster( nCores ) 
  registerDoParallel(cl.2)
  clusterExport(cl.2, "raw.data.longitude")
  interpSurfaceMetric[ , "Lon" ] <- parApply(cl.2, interpSurfaceMetric[,c("Var1","Var2")], 1 , function(x) raw.data.longitude[x[1],x[2]])
  stopCluster(cl.2) ; rm(cl.2); gc(reset=TRUE)
  
  cl.2 <- makeCluster( nCores ) 
  registerDoParallel(cl.2)
  clusterExport(cl.2, "raw.data.latitude")
  interpSurfaceMetric[ , "Lat" ] <- parApply(cl.2, interpSurfaceMetric[,c("Var1","Var2")], 1 , function(x) raw.data.latitude[x[1],x[2]])
  stopCluster(cl.2) ; rm(cl.2); gc(reset=TRUE)
  
}

if( ! multdimCoordinates ) {
  
  cl.2 <- makeCluster( nCores ) 
  registerDoParallel(cl.2)
  clusterExport(cl.2, "raw.data.longitude")
  interpSurfaceMetric[ , "Lon" ] <- parSapply(cl.2, interpSurfaceMetric[,"Var1"] , function(x) raw.data.longitude[x])
  stopCluster(cl.2) ; rm(cl.2); gc(reset=TRUE)
  
  cl.2 <- makeCluster( nCores ) 
  registerDoParallel(cl.2)
  clusterExport(cl.2, "raw.data.latitude")
  interpSurfaceMetric[ , "Lat" ] <- parSapply(cl.2, interpSurfaceMetric[,"Var2"] , function(x) raw.data.latitude[x])
  stopCluster(cl.2) ; rm(cl.2); gc(reset=TRUE)
  
}

if( rawDataIs3D ) { 
  
  cl.2 <- makeCluster( nCores ) 
  registerDoParallel(cl.2)
  clusterExport(cl.2, "raw.data.depth")
  interpSurfaceMetric[ , "Depth" ] <- parSapply(cl.2, interpSurfaceMetric[,"Var3"] , function(x) raw.data.depth[x])
  stopCluster(cl.2) ; rm(cl.2); gc(reset=TRUE)
  
}

## --------------------------------------------------------

shape <- raster(paste0(bathymetryFolder,"/",bathymetryFiles[1]))
shape <- crop(shape,extent(outcomeExtent))
shape[!is.na(shape)] <- 1

rawDataCoordsTest <- interpSurfaceMetric[interpSurfaceMetric$Var3 == 1 , c("Lon","Lat") ]
rawDataResolution <- max( c(360/dim.i,180/dim.j) )
rawDataTest <- rasterFromXYZ.V2(rawDataCoordsTest,rawDataVarTest)

## ----------------

png(filename = paste0(processTempFolder,"/rawDataTestImage.png"), width = 1920, height = 1080)
plot(rawDataTest)
dev.off()

writeRaster(rawDataTest,filename=paste0(processTempFolder,"/rawDataTestImage.tif"),format="GTiff",overwrite=T)

## --------------------------------------------------------
## --------------------------------------------------------

file.remove( list.files(tempFolder, full.names = TRUE, pattern = "depthSteps") )
file.remove( list.files(tempFolder, full.names = TRUE, pattern = "rawData") ) 

# --------------------------

if( ! is.na(dim(raw.data.longitude)[2]) ) { raw.data.longitudeDim <- 1:dim(raw.data.longitude)[1] } else { raw.data.longitudeDim <- 1:length(raw.data.longitude) }
if( ! is.na(dim(raw.data.latitude)[2]) ) { raw.data.latitudeDim <- 1:dim(raw.data.latitude)[2] } else { raw.data.latitudeDim <- 1:length(raw.data.latitude) }

timeListAvailableYDaysSince <- as.numeric(as.Date(paste0(timeListAvailableY,"-01-01"), "%Y-%m-%d"))

ncDim.Lon <- ncdim_def( "Lon", "degreesE", raw.data.longitudeDim)
ncDim.Lat <- ncdim_def( "Lat", "degreesN", raw.data.latitudeDim)
ncDim.Depth <- ncdim_def( "Depth", "meter", 1:2)
ncDim.Time <- ncdim_def( "Time", "days since 1970-01-01", timeListAvailableYDaysSince, unlim=TRUE)
raw.dataMaxVar <- ncvar_def("raw.dataMax", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
raw.dataMeanVar <- ncvar_def("raw.dataMean", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
raw.dataMinVar <- ncvar_def("raw.dataMin", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
ncDims <- list(raw.dataMinVar,raw.dataMeanVar,raw.dataMaxVar)

# --------------------------

seqListing <- c(raw.data.depth)
seqChunks <- data.frame(from = c(0,seqListing[-length(seqListing)]), to = c(seqListing[-length(seqListing)],11000) )

if( exists("truncateRawDataToDepthLayers") ) { if( truncateRawDataToDepthLayers ) { 

  interpSurfaceMetric <- interpSurfaceMetric[ which(interpSurfaceMetric$Depth %in% raw.data.depth[sapply(depthsToInterpolate, function(x) { which.min(abs(seqChunks[,1] - x)) } )]) ,]
  seqChunks <- seqChunks[sapply(depthsToInterpolate, function(x) { which.min(abs(seqChunks[,1] - x)) } ),]

}}

parallelSequence <- expand.grid(depthRange=1:nrow(seqChunks))

if( rawDataIs3D & realm.i == "Surface" ) { parallelSequence <- data.frame(depthRange=parallelSequence[parallelSequence$depthRange == 1 , ]) }
if( ! rawDataIs3D ) { parallelSequence <- seqChunks }

## ------------------------------

Cluster <- makeCluster(nCores)
registerDoParallel(Cluster) 

producefiles <- foreach(d=1:nrow(parallelSequence), .verbose=FALSE, .noexport = c("arma3DMin","arma3DMean","arma3DMax"),  .packages=c("stringr","Rcpp","ncdf4.helpers","abind","ncdf4","data.table","bigmemory","easyNCDF")) %dopar% {
  
  # ----------------------
  
  source("Dependencies/mainFunctionsArmadillo.R")
  
  # ----------------------
  
  nc.file <- paste0(tempFolder,"/","rawDataDepth_",d,".nc")
  if( file.exists(nc.file) ) { file.remove(nc.file) }
  nc <- nc_create( nc.file, ncDims )
  nc_close(nc)
  
  # ----------------------
  
  dQ <- which(raw.data.depth %in% unlist(seqChunks[d,]) )
  
  if( d == 1 ) { dQ <- c(1,1) }
  if( d == nrow(seqChunks) & d != 1 ) { dQ <- c(nrow(seqChunks)-1,nrow(seqChunks)) }
  if( rawDataIs3D ) { dQ.i.loop <- 1:2 } else { dQ.i.loop <- 1 }
  
  # ----------------------
  
  for( y in timeListAvailableY ) {
    
    for( dQ.i in dQ.i.loop ) {
        
        ## --------------------
      
        temArrayDepth <- array(dim = c(dim.i, dim.j, 12))
        noDataVector <- numeric(0)
        
        for( m in 1:12) {

          nc.file <- timeListAvailable[ which(substr(timeListAvailable$period,1,4) %in% y & as.numeric(substr(timeListAvailable$period,6,7)) %in% m) , "file"]
          
          if(length(nc.file) == 0) { noDataVector <- c(noDataVector,m); next }

          raw.data <- nc_open( nc.file, readunlim=FALSE )
          
          errorCatch <- FALSE
          tryCatch( nc.fileTime <- nc.get.time.series(raw.data), error = function(e){ errorCatch <<- TRUE } )    
          if( ! errorCatch ) {  nc.fileTime <- as.Date(as.character(nc.fileTime)) }
          if( errorCatch ) { nc.fileTime <- as.Date(ncvar_get( raw.data, "time"),origin="1900-01-01") }

          if(is.na(nc.fileTime) & singleFile ) { 
            nc.fileTime <- substr(nc.file , unlist(gregexpr(variable.i,nc.file)) + nchar(variable.i) + 1 , unlist(gregexpr("-",nc.file)) - 1 ) 
            nc.fileTime <- as.Date(paste0(substr(nc.fileTime,1,4),"-",str_pad(substr(nc.fileTime,5,6), 2, pad = "0"),"-",str_pad(substr(nc.fileTime,8,9), 2, pad = "0")))
          }
          
          loc.t <- which( as.numeric(substr(nc.fileTime,1,4)) == y & as.numeric(substr(nc.fileTime,6,7)) == m )
          
          if(   rawDataIs3D ) { temArrayDepth[,,m] <- ncvar_get( raw.data, variable.i , start=c(1,1,dQ[dQ.i],loc.t), count=c(dim.i,dim.j,1,1) ) } 
          if( ! rawDataIs3D & ! singleFile ) { temArrayDepth[,,m] <- ncvar_get( raw.data, variable.i , start=c(1,1,loc.t), count=c(dim.i,dim.j,1) ) } 
          if( ! rawDataIs3D & singleFile) { temArrayDepth[,,m] <- ncvar_get( raw.data, variable.i , start=c(1,1), count=c(dim.i,dim.j) ) } 
          
          if( ! is.null(correct.raw.data.matrix) ) {
            
              for( ii in 1:nrow(correct.raw.data.matrix)) {
                for( jj in 1:ncol(correct.raw.data.matrix)) {
                  if( ! is.na(correct.raw.data.matrix[ii,jj]) & is.na( temArrayDepth[ii,jj,m] ) ) { temArrayDepth[ii,jj,m] <- 0 }
                }
              }
            
          }
          
          nc_close( raw.data )
          
        }
          
        if( length(noDataVector) == 12 ) { stop("Error :: No data collected from nc files")}
        if( length(noDataVector) > 0 ) { temArrayDepth <- temArrayDepth[,,(1:12)[-which(1:12 %in% noDataVector)]]}
        
        ## --------------------

        temArrayMeanDepth.i <- arma3DMean(temArrayDepth,3)
        temArrayMinDepth.i <- arma3DMin(temArrayDepth,3)
        temArrayMaxDepth.i <- arma3DMax(temArrayDepth,3)

        ## --------------------
        
        if(dQ.i == 1) {
          
          temArrayMeanDepth <- temArrayMeanDepth.i
          temArrayMinDepth <- temArrayMinDepth.i
          temArrayMaxDepth <- temArrayMaxDepth.i

    
        }
        
        ## --------------------
        
        if(dQ.i == 2) {

          temArrayMeanDepth <- abind( temArrayMeanDepth, temArrayMeanDepth.i, along=3 )
          temArrayMinDepth <- abind( temArrayMinDepth, temArrayMinDepth.i, along=3 )
          temArrayMaxDepth <- abind( temArrayMaxDepth, temArrayMaxDepth.i, along=3 )
          
        }
      
        ## --------------------
        
        rm(temArrayDepth)

    }
    
    ## --------------------
    
    if( ! rawDataIs3D ) { 
    
      temArrayMeanDepth <- abind(temArrayMeanDepth,temArrayMeanDepth,along=3)
      temArrayMinDepth <- abind(temArrayMinDepth,temArrayMinDepth,along=3)
      temArrayMaxDepth <- abind(temArrayMaxDepth,temArrayMaxDepth,along=3)
      
    }
    
    ## --------------------

    t.loc <- which(timeListAvailableY == y)

    nc.file <- paste0(tempFolder,"/","rawDataDepth_",d,".nc")
    nc <- nc_open( nc.file, readunlim=FALSE , write=TRUE )
    ncvar_put(nc, varid = "raw.dataMin", vals = temArrayMinDepth, start = c(1, 1, 1, t.loc), count = c(dim(temArrayMinDepth)[1], dim(temArrayMinDepth)[2], dim(temArrayMinDepth)[3],1) )
    ncvar_put(nc, varid = "raw.dataMean", vals = temArrayMeanDepth, start = c(1, 1, 1, t.loc), count = c(dim(temArrayMeanDepth)[1], dim(temArrayMeanDepth)[2], dim(temArrayMeanDepth)[3],1) )
    ncvar_put(nc, varid = "raw.dataMax", vals = temArrayMaxDepth, start = c(1, 1, 1, t.loc), count = c(dim(temArrayMaxDepth)[1], dim(temArrayMaxDepth)[2], dim(temArrayMaxDepth)[3],1) )
    nc_close( nc )
    
    rm(temArrayMeanDepth); rm(temArrayMinDepth); rm(temArrayMaxDepth)
    gc(reset=TRUE)
    
  }
  
  ## --------------------
  ## --------------------
  ## If to be anomaly
  
  if( experiment != "Baseline" ) {
    
    timeListAvailableDaysSince <- as.numeric(as.Date(paste0(outcomePeriodY.available,"-01-01"), "%Y-%m-%d"))
    ncDim.Time <- ncdim_def( "Time", "days since 1970-01-01", timeListAvailableDaysSince, unlim=TRUE)
    raw.dataMaxVar <- ncvar_def("raw.dataMax", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
    raw.dataMeanVar <- ncvar_def("raw.dataMean", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
    raw.dataMinVar <- ncvar_def("raw.dataMin", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
    ncDims <- list(raw.dataMinVar,raw.dataMeanVar,raw.dataMaxVar)
    
    nc.file <- paste0(tempFolder,"/","rawDataDepth_",d,".nc")
    nc <- nc_open( nc.file, readunlim=FALSE , write=TRUE )
    
    nc.file.dec <- paste0(tempFolder,"/","rawDataTempDepth_",d,".nc")
    if( file.exists(nc.file.dec) ) { file.remove(nc.file.dec) }
    ncNew <- nc_create( nc.file.dec, ncDims )
    
    t.loc <- which(timeListAvailableY %in% baselinePeriod)
    
    tempArrayBaselineMin <- ncvar_get( nc, "raw.dataMin", start=c(1,1,1,t.loc[1]), count=c(dim.i,dim.j,2,t.loc[length(t.loc)] - t.loc[1] + 1))
    tempArrayBaselineMin <- means.along(tempArrayBaselineMin, 4)
    
    tempArrayBaselineMean <- ncvar_get( nc, "raw.dataMean", start=c(1,1,1,t.loc[1]), count=c(dim.i,dim.j,2,t.loc[length(t.loc)] - t.loc[1] + 1))
    tempArrayBaselineMean <- means.along(tempArrayBaselineMean, 4)
    
    tempArrayBaselineMax <- ncvar_get( nc, "raw.dataMax", start=c(1,1,1,t.loc[1]), count=c(dim.i,dim.j,2,t.loc[length(t.loc)] - t.loc[1] + 1))
    tempArrayBaselineMax <- means.along(tempArrayBaselineMax, 4)
    
    t.loc.i <- which(timeListAvailableY %in% outcomePeriodY)
    
    for( t.i in t.loc.i) {
      
      for(d.i in 1:2) {
    
        t.loc.i.new <- which(outcomePeriodY.available == timeListAvailableY[t.i])    
        tempArrayProjection <- ncvar_get( nc, "raw.dataMin", start=c(1,1,d.i,t.i), count=c(dim.i,dim.j,1,1))
        tempArrayProjection <- sweep(tempArrayProjection,1,tempArrayBaselineMin[,,d.i], FUN = "-")
        ncvar_put(ncNew, varid = "raw.dataMin", vals = tempArrayProjection, start = c(1, 1, d.i, t.loc.i.new), count = c(dim.i, dim.j, 1, 1) )
        
        tempArrayProjection <- ncvar_get( nc, "raw.dataMean", start=c(1,1,d.i,t.i), count=c(dim.i,dim.j,1,1))
        tempArrayProjection <- sweep(tempArrayProjection,1,tempArrayBaselineMean[,,d.i], FUN = "-")
        ncvar_put(ncNew, varid = "raw.dataMean", vals = tempArrayProjection, start = c(1, 1, d.i, t.loc.i.new), count = c(dim.i, dim.j, 1, 1) )
        
        tempArrayProjection <- ncvar_get( nc, "raw.dataMax", start=c(1,1,d.i,t.i), count=c(dim.i,dim.j,1,1))
        tempArrayProjection <- sweep(tempArrayProjection,1,tempArrayBaselineMax[,,d.i], FUN = "-")
        ncvar_put(ncNew, varid = "raw.dataMax", vals = tempArrayProjection, start = c(1, 1, d.i, t.loc.i.new), count = c(dim.i, dim.j, 1, 1) )
        
      }
    }

    nc_close(nc)
    nc_close(ncNew)
    
    file.remove(nc.file)
    file.rename(nc.file.dec,nc.file)
    
  }
  
  ## --------------------
  ## --------------------
  ## If to be decadal
  
  if( outcomePeriodType == "decade" ) {
    
    ncDim.Time <- ncdim_def( "Time", "days since 1970-01-01", timeListAvailableDecDaysSince, unlim=TRUE)
    raw.dataMaxVar <- ncvar_def("raw.dataMax", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
    raw.dataLtMaxVar <- ncvar_def("raw.dataLtMax", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
    raw.dataMeanVar <- ncvar_def("raw.dataMean", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
    raw.dataLtMinVar <- ncvar_def("raw.dataLtMin", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
    raw.dataMinVar <- ncvar_def("raw.dataMin", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
    raw.dataRangeVar <- ncvar_def("raw.dataRange", "unit",  list(ncDim.Lon,ncDim.Lat,ncDim.Depth,ncDim.Time), 1.e30 ) 
    ncDims <- list(raw.dataMinVar,raw.dataMeanVar,raw.dataMaxVar,raw.dataLtMaxVar,raw.dataLtMinVar,raw.dataRangeVar)
    
    nc.file <- paste0(tempFolder,"/","rawDataDepth_",d,".nc")
    nc <- nc_open( nc.file, readunlim=FALSE , write=TRUE )
    
    nc.file.dec <- paste0(tempFolder,"/","rawDataDecDepth_",d,".nc")
    if( file.exists(nc.file.dec) ) { file.remove(nc.file.dec) }
    ncNew <- nc_create( nc.file.dec, ncDims )
    
    for( dec.i in 1:nrow(decadalTimeSteps) ) {
      
      t.loc <- which(outcomePeriodY.available %in% decadalTimeSteps[dec.i,1]:decadalTimeSteps[dec.i,2])
      
      for(d.i in 1:2) {
        
        temArrayToDecade <- ncvar_get( nc, "raw.dataMax", start=c(1,1,d.i,t.loc[1]), count=c(dim.i,dim.j,1,t.loc[length(t.loc)] - t.loc[1] + 1))
        ncvar_put(ncNew, varid = "raw.dataMax", vals = arma3DMax(temArrayToDecade,3) , start = c(1, 1, d.i, dec.i), count = c(dim(temArrayToDecade)[1], dim(temArrayToDecade)[2], 1 , 1) )
        ncvar_put(ncNew, varid = "raw.dataLtMax", vals = arma3DMean(temArrayToDecade,3) , start = c(1, 1, d.i, dec.i), count = c(dim(temArrayToDecade)[1], dim(temArrayToDecade)[2], 1 , 1) )
        
        temArrayDecadeMax <- arma3DMax(temArrayToDecade,3)
        
        temArrayToDecade <- ncvar_get( nc, "raw.dataMean", start=c(1,1,d.i,t.loc[1]), count=c(dim.i,dim.j,1,t.loc[length(t.loc)] - t.loc[1] + 1))
        ncvar_put(ncNew, varid = "raw.dataMean", vals = arma3DMean(temArrayToDecade,3) , start = c(1, 1, d.i, dec.i), count = c(dim(temArrayToDecade)[1], dim(temArrayToDecade)[2], 1 , 1) )
        
        temArrayToDecade <- ncvar_get( nc, "raw.dataMin", start=c(1,1,d.i,t.loc[1]), count=c(dim.i,dim.j,1,t.loc[length(t.loc)] - t.loc[1] + 1))
        ncvar_put(ncNew, varid = "raw.dataMin", vals = arma3DMin(temArrayToDecade,3) , start = c(1, 1, d.i, dec.i), count = c(dim(temArrayToDecade)[1], dim(temArrayToDecade)[2], 1 , 1) )
        ncvar_put(ncNew, varid = "raw.dataLtMin", vals = arma3DMean(temArrayToDecade,3) , start = c(1, 1, d.i, dec.i), count = c(dim(temArrayToDecade)[1], dim(temArrayToDecade)[2], 1 , 1) )
        
        temArrayDecadeMin <- arma3DMin(temArrayToDecade,3)
        
        ncvar_put(ncNew, varid = "raw.dataRange", vals = temArrayDecadeMax - temArrayDecadeMin , start = c(1, 1, d.i, dec.i), count = c(dim(temArrayToDecade)[1], dim(temArrayToDecade)[2], 1 , 1) )
        
      }
    }
    
    nc_close(ncNew)
    nc_close(nc)
    
    file.remove(nc.file)
    file.rename(nc.file.dec,nc.file)
    
  }
  
  ## --------------------
  ## --------------------
  
  return(NULL)

}

stopCluster(Cluster) ; rm(Cluster)

## --------------------------------------------------------
## --------------------------------------------------------
# Subset by final product extent

interpSurfaceMetric <- interpSurfaceMetric[which(interpSurfaceMetric$Lon >= (outcomeExtent[1] - 5) & interpSurfaceMetric$Lon <= (outcomeExtent[2] + 5) & interpSurfaceMetric$Lat >= (outcomeExtent[3] - 5) & interpSurfaceMetric$Lat <= (outcomeExtent[4] + 5)) ,]

## --------------------------------------------------------
## --------------------------------------------------------
# Subset by final product extent

metadata <- data.frame(   descriptor = c("variableLongName" ,
                                         "variableStandardName",
                                         "variable",
                                         "rawDatais3D",
                                         "rawDataResolution",
                                         "units",
                                         "dimX",
                                         "dimY",
                                         "dimZ",
                                         "climatologyType",
                                         "dimT",
                                         "minT",
                                         "maxT",
                                         "dataSource",
                                         "dataProduct",
                                         "dataDataset",
                                         "minValues",
                                         "maxValues"),
                          value = c(variableLongName,
                                    variableStandardName,
                                    variable,
                                    rawDataIs3D,
                                    paste(c(rawDataResolution,dim.z), collapse=","),
                                    metadataUnits,
                                    dim.i,
                                    dim.j,
                                    ifelse(is.null(dim.z),"NULL",dim.z),
                                    outcomePeriodType,
                                    dim.t,
                                    min(outcomePeriodY.available),
                                    max(outcomePeriodY.available),
                                    ifelse(is.null(metadataSource),"NA",metadataSource),
                                    ifelse(is.null(metadataProduct),"NA",metadataProduct),
                                    metadataDataset,
                                    cellStats(rawDataTest,min,na.rm=T),
                                    cellStats(rawDataTest,max,na.rm=T)) )

metadata$value <- gsub(";", ",", metadata$value, fixed = TRUE)
write.table(metadata,file=paste0(processTempFolder,"/rawDataMetadata.txt"),col.names = TRUE,sep=";",dec=".",row.names = FALSE,quote=FALSE)

## --------------------------------------------------------
## --------------------------------------------------------

rawDataReady <- TRUE
rm( list = c("shape","rawDataVarTest","rawDataTest","raw.data.longitude","raw.data.latitude","raw.data.var","raw.data"))

closeAllConnections()
gc(reset=TRUE, full = TRUE)

## ------------------ 

cat("# End process // ET:", round(difftime(Sys.time(), time.i, units='hours'), digits = 3) ,"hours\n")

## --------------------------------------------------------
## --------------------------------------------------------