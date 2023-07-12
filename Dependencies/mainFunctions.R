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

## ----------------------
## ----------------------

packagesToUse <- c("rWind","BiocManager","Rcpp","bigmemory","ggplot2","abind","nctools","RNetCDF","ff","thredds","reshape2","jsonlite","matrixStats","easyNCDF","FNN","raster","ncdf4","data.table","matrixStats","gstat","bigmemory","doParallel","parallel","compiler","ncdf4.helpers")

for(package in packagesToUse) {
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package)}
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package , type = "source") }
  if( ! package %in% rownames(installed.packages()) ) { BiocManager::install(package) }
  if( ! package %in% rownames(installed.packages()) & package == "thredds" ) { devtools::install_github("bocinsky/thredds") }
  if( ! package %in% rownames(installed.packages()) & package == "nctools" ) { devtools::install_github("roliveros-ramos/nctools") }
  
  if( ! package %in% rownames(installed.packages()) ) { stop("Error on package instalation") }
  library(package, character.only = TRUE)
}

## ---------------------------
## ---------------------------

source("Dependencies/mainFunctionsArmadillo.R")

## ------------------------

log10Ceiling <- function(x) { 10^(ceiling(log10(x))) }

## ---------------------------

firstToupper <- function(string) {
  
  return(paste(toupper(substr(string, 1, 1)), substr(string, 2, nchar(string)), sep=""))

}

## ---------------------------

interpNARaster <- function(raster,shape) {
      
    dataCells <- as.data.frame(raster, xy=TRUE, na.rm=TRUE)
    interpCells <- Which(is.na(raster), cells=TRUE)[!Which(is.na(raster), cells=TRUE) %in% Which(is.na(shape), cells=TRUE)]
    
    interpCellsLoc <- data.frame(x=sapply(interpCells,function(x) xyFromCell(raster,x)[,1]),
                                 y=sapply(interpCells,function(x) xyFromCell(raster,x)[,2]))
    
    idw.nearest.r <- get.knnx( dataCells[ , c("x","y") ] , as.matrix(interpCellsLoc) , k = 3 , algorithm="kd_tree" )
    idw.nearest.i <- idw.nearest.r$nn.index
    idw.nearest.d <- idw.nearest.r$nn.dist
    interpCellsValue <- numeric(length(interpCells))
    
    for( c in 1:length(interpCells)) {
      
      interpCellsValue[c] <- (sum( dataCells[idw.nearest.i[c,],3] / idw.nearest.d[c,]^idwPower , na.rm=TRUE)) / (sum( 1 / idw.nearest.d[c,]^idwPower , na.rm=TRUE))
      
    }
    
    raster[interpCells] <- interpCellsValue
    return(raster)
}

## ---------------------------

means.along <- function(a, i) {
  n <- length(dim(a))
  b <- aperm(a, c(seq_len(n)[-i], i))
  rowMeans(b, dims = n - 1)
}

## ---------------------------

listMemory <- cmpfun( function (pos = 1, pattern, order.by = "Size", decreasing=TRUE, head = TRUE, n = 100) {
  
  napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = pos)))
  
  names <- ls(pos = pos, pattern = pattern)
  
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  
  obj.mode <- napply(names, mode)
  
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  
  obj.size <- napply(names, object.size) / 10^6 # megabytes
  
  obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
  
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  
  obj.dim[vec, 1] <- napply(names, length)[vec]
  
  out <- data.frame(obj.type, obj.size, obj.dim)
  
  names(out) <- c("Type", "Size", "Rows", "Columns")
  
  out <- out[order(out[[order.by]], decreasing=decreasing), ]
  
  if (head)
    
    out <- head(out, n)
  
  out
  
  
  
} )

## ---------------------------

ListingAvailableRawData <- function(folder) {
  
  filesAvailableName <- list.files(folder, full.names=F, pattern=".nc",recursive=T)
  filesAvailableURL <- list.files(folder, full.names=T, pattern=".nc",recursive=T)
  
  availableRawData <- data.frame()
  
  for( i in 1:length(filesAvailableName)) {
    
    filesAvailableName.i <- filesAvailableName[i]
    fileBreaks <- unlist(gregexpr("_",filesAvailableName.i))
    
    variable <- substr(filesAvailableName.i,1,fileBreaks[1]-1)
    frequency <- substr(filesAvailableName.i,fileBreaks[1]+1,fileBreaks[2]-1)
    model <- substr(filesAvailableName.i,fileBreaks[2]+1,fileBreaks[3]-1)
    experiment <- substr(filesAvailableName.i,fileBreaks[3]+1,fileBreaks[4]-1)
    ensemble <- substr(filesAvailableName.i,fileBreaks[4]+1,fileBreaks[5]-1)
    
    fileBreaks <- unlist(gregexpr("-",filesAvailableName.i))
    fileBreaks <- fileBreaks[length(fileBreaks)]
    period <- substr(filesAvailableName.i,fileBreaks-6,fileBreaks+6)

    if( exists("dimNameLon") ) { rm(dimNameLon, pos = ".GlobalEnv"); rm(dimNameLat, pos = ".GlobalEnv"); rm(dimNameDepth, pos = ".GlobalEnv") }
    
    namegeo <- "NA"
    nameDepth <- "NA"
    validFile <- TRUE
    
    tryCatch( getDimNames(filesAvailableURL[i]) , error=function(e){ validFile <<- FALSE } )
    
    if( validFile ) {
      
      namegeo <- getDimUnits(filesAvailableURL[i],dimNameLon)
      tryCatch( nameDepth <- getDimUnits(filesAvailableURL[i],dimNameDepth) , error=function(e){ nameDepth <<- "NA" } )
      
    }

    availableRawData.i <- data.frame(variable=variable,
                                     frequency=frequency,
                                     model=model,
                                     experiment=experiment,
                                     ensemble=ensemble,
                                     period=period,
                                     dimNameGeo=namegeo,
                                     dimNameDepth=nameDepth,
                                     validFile=validFile
                                     )
    
    if( file.info(filesAvailableURL[i])$size > 100000) { availableRawData <- rbind(availableRawData,availableRawData.i) }
    
  }

  availableRawData <- availableRawData[ which(availableRawData$validFile), ]
  availableRawData <- unique(availableRawData)
  
  availableRawDataOverall <- availableRawData[!duplicated(availableRawData[,1:5]),]
  availableRawDataOverall$periodFrom <- NA
  availableRawDataOverall$periodTo <- NA

  for( i in 1:nrow(availableRawDataOverall)) {
    
    machingFiles <- which(availableRawData$variable == availableRawDataOverall[i,"variable"] &
          availableRawData$frequency == availableRawDataOverall[i,"frequency"] &
          availableRawData$model == availableRawDataOverall[i,"model"] &
          availableRawData$experiment == availableRawDataOverall[i,"experiment"] &
          availableRawData$ensemble == availableRawDataOverall[i,"ensemble"] )
          
    availableRawDataOverall.i <- unique(sort(unlist(apply(data.frame(from = as.numeric(substr(availableRawData[machingFiles,"period"],1,4)),to = as.numeric(substr(availableRawData[machingFiles,"period"],8,11))), 1 , function(x) { x[1]:x[2] } ))))
    
    if( class(availableRawDataOverall.i) != "integer") { stop("!")}
    
    if(length(min(availableRawDataOverall.i):max(availableRawDataOverall.i)) == length(availableRawDataOverall.i)) { 
      
      availableRawDataOverall[i,"periodFrom"] <- min(availableRawDataOverall.i)
      availableRawDataOverall[i,"periodTo"] <- max(availableRawDataOverall.i)
      
      }
  
    dimNameDepth <- unique(availableRawData[machingFiles,"dimNameDepth"])
    dimNameDepth <- ifelse( length(dimNameDepth) > 1 , dimNameDepth[dimNameDepth != "NA"] , dimNameDepth )
    availableRawDataOverall[i,"dimNameDepth"] <- dimNameDepth[1]
    
  }
  
  availableRawDataOverall <- availableRawDataOverall[which(!is.na(availableRawDataOverall$periodFrom)),]
  availableRawDataOverall <-  availableRawDataOverall[with(availableRawDataOverall, order(variable, experiment,periodFrom)),]
  
  return(availableRawDataOverall)

}

## ---------------------------

ListingZeroSizeRawData <- function(folder) {
  
  filesAvailableName <- list.files(folder, full.names=F, pattern=".nc",recursive=T)
  filesAvailableURL <- list.files(folder, full.names=T, pattern=".nc",recursive=T)
  
  availableRawData <- data.frame()
  
  for( i in 1:length(filesAvailableName)) {
    
    fileSize <- file.info(filesAvailableURL[i])$size
    
    filesAvailableName.i <- filesAvailableName[i]
    fileBreaks <- unlist(gregexpr("_",filesAvailableName.i))
    
    variable <- substr(filesAvailableName.i,1,fileBreaks[1]-1)
    frequency <- substr(filesAvailableName.i,fileBreaks[1]+1,fileBreaks[2]-1)
    model <- substr(filesAvailableName.i,fileBreaks[2]+1,fileBreaks[3]-1)
    experiment <- substr(filesAvailableName.i,fileBreaks[3]+1,fileBreaks[4]-1)
    ensemble <- substr(filesAvailableName.i,fileBreaks[4]+1,fileBreaks[5]-1)
    
    fileBreaks <- unlist(gregexpr("-",filesAvailableName.i))
    fileBreaks <- fileBreaks[length(fileBreaks)]
    period <- substr(filesAvailableName.i,fileBreaks-6,fileBreaks+6)
    
    availableRawData.i <- data.frame(variable=variable,
                                     frequency=frequency,
                                     model=model,
                                     experiment=experiment,
                                     ensemble=ensemble,
                                     period=period,
                                     fileSize=fileSize
    )
    
    availableRawData <- rbind(availableRawData,availableRawData.i)
    
  }
  
  if( sum(availableRawData$fileSize == 0) > 0 ) {
    
    cat(sum(availableRawData$fileSize == 0),"file(s) with 0 size","\n")
    availableRawDataOverall <- availableRawData[which(availableRawData$fileSize == 0),]
    availableRawDataOverall <- availableRawDataOverall[!duplicated(availableRawDataOverall[,1:5]),]
    return(availableRawDataOverall)
    
  }
  if( sum(availableRawData$fileSize == 0) == 0 ) {
    
    cat("No file with 0 size \n")
    return(NULL)
    
  }
 
}

## ---------------------------

ListingCMIP <- function(model,experiment,realm,variable,timeFreq) {
  
  file.remove(list.files("." , pattern=paste0( "query.json") , full.names=TRUE))
  jsonData <- numeric(0)
  
  try(download.file(paste0("https://esgf-node.llnl.gov/esg-search/search?latest=true&variable=",variable,"&experiment_id=",experiment,"&limit=10000&format=application%2Fsolr%2Bjson"),destfile=paste0("query.json"),method="auto",quiet=TRUE), silent = TRUE)
  
  jsonFile <- readLines( paste0( "query.json") )
  
  if( file.exists(paste0( "query.json")) & ! length(jsonFile) == 0 ) {
    
    jsonData <- jsonlite::fromJSON( jsonFile )
    jsonData <- jsonData$response
    jsonData <- jsonData$docs
    
  }
  
  file.remove(list.files("." , pattern=paste0( "query.json") , full.names=TRUE))
  
  if( length(jsonData) > 0) { 
    
    jsonData <- jsonData[which(unlist(jsonData$source_id) %in% model),]
    jsonData <- jsonData[which(unlist(jsonData$table_id) %in% timeFreq),]
    jsonData <- jsonData[which(unlist(jsonData$realm) %in% realm),]
    
  }
  
  return(jsonData)
}

## ---------------------------

downloadCMIP <- function(listProject,folder,period,concurrent,overwriteFiles) {
  
  variable <- listProject$variable[[1]]
  
  comb <- data.frame( source_id=unlist(listProject$source_id) , member_id=unlist(listProject$member_id))
  comb <- data.frame(comb,variable=unique(unlist(listProject$variable)),experiment_id=unique(unlist(listProject$experiment_id)),table_id=unique(unlist(listProject$table_id)))
  comb <- unique(comb)
  
  finalListModels <- unique(unlist(comb$source_id))
  
  cl.2 <- parallel::makeCluster(concurrent)
  registerDoParallel(cl.2)
  
  downloads <- foreach(m = 1:length(finalListModels), .verbose=FALSE, .packages=c("bigmemory")) %dopar% {
    
    comb.i <- comb[unlist(comb$source_id) == finalListModels[m],]
    finalMember_id <- unique(unlist(comb.i$member_id))
    finalMember_id <- finalMember_id[sort(as.numeric(as.factor(finalMember_id)),index.return=T)$ix]
    finalMember_id <- c(finalMember_id[which(grepl("r1i",finalMember_id))],finalMember_id)[1]
    
    for( f in 1:length(finalMember_id)) {
      
      # If to overwrite
      if( overwriteFiles & length(list.files(folder,pattern=paste0(variable,"_",comb.i[1,"table_id"],"_",comb.i[1,"source_id"],"_",comb.i[1,"experiment_id"]))) > 0 ) { file.remove(list.files(folder,pattern=paste0(variable,"_",comb.i[1,"table_id"],"_",comb.i[1,"source_id"],"_",comb.i[1,"experiment_id"]), full.names = TRUE)) }
      
      # Remove those with no information
      filesToRemove <- which(sapply(list.files(folder,pattern=paste0(variable,"_",comb.i[1,"table_id"],"_",comb.i[1,"source_id"],"_",comb.i[1,"experiment_id"]), full.names = TRUE), function(x) { file.info(x)$size } ) < 1000)
      if( length(filesToRemove) > 0 ) { file.remove( list.files(folder,pattern=paste0(variable,"_",comb.i[1,"table_id"],"_",comb.i[1,"source_id"],"_",comb.i[1,"experiment_id"]), full.names = TRUE)[filesToRemove]) }
      
      # Full listed and no need to download
      if( ! overwriteFiles & length(list.files(folder,pattern=paste0(variable,"_",comb.i[1,"table_id"],"_",comb.i[1,"source_id"],"_",comb.i[1,"experiment_id"]))) > 0 ) { next }
      
      # -------------
      
      fileURL <- paste0("https://esgf-node.llnl.gov/esg-search/wget?variable=",comb.i[1,"variable"],"&experiment_id=",comb.i[1,"experiment_id"],"&source_id=",comb.i[1,"source_id"],"&table_id=",comb.i[1,"table_id"],"&member_id=",finalMember_id[f])
      
      passError <- TRUE
      attempts <- 0
      
      while( passError ) {
        passError <- FALSE
        attempts <- attempts + 1
        tryCatch( download.file(fileURL,destfile=paste0(m,"_",f,"wget.sh"), method ="curl",quiet=TRUE, timeout=100000000,extra='-L -k') , error=function(e){ passError <<- TRUE } )
        if( passError ) { tryCatch( download.file(fileURL,destfile=paste0(m,"_",f,"wget.sh")), error=function(e){ passError <<- TRUE } ) }
        if(attempts > 9) { break }
      }

      if( passError ) { next }
      
      wgetFile <- readLines(paste0(m,"_",f,"wget.sh"))
      file.remove(list.files("." , pattern=paste0(m,"_",f,"wget.sh") , full.names=TRUE))
      
      lines.in.parse <- which(sapply(1:length(wgetFile),function(x) { grepl(".nc",wgetFile[x]) & grepl("http://",wgetFile[x]) }))
      
      if( length(lines.in.parse) == 0) { lines.in.parse <- which(sapply(1:length(wgetFile),function(x) { grepl(".nc",wgetFile[x]) & grepl("https://",wgetFile[x]) })) }
      
      if( length(lines.in.parse) == 0) { next }
      
      for( f.i in (1:length(lines.in.parse))[(1:length(lines.in.parse) != 0)] ) {
        
        wgetFile.i <- wgetFile[lines.in.parse[f.i]]
        
        if( grepl("LGM",toupper(gsub("'","",strsplit(wgetFile.i," ")[[1]][1]))) | grepl("LIG",toupper(gsub("'","",strsplit(wgetFile.i," ")[[1]][1]))) | grepl("MH",toupper(gsub("'","",strsplit(wgetFile.i," ")[[1]][1]))) ) { stop("Error :: 510") }
        
        # Test for data within period
        testForPeriod <- gsub("'","",strsplit(wgetFile.i," ")[[1]][1])
        
        if(nchar(testForPeriod) == 1 ) { next }
        
        testForPeriodLoc.1 <- unlist(gregexpr("_",testForPeriod))[length(unlist(gregexpr("_",testForPeriod)))]
        testForPeriodLoc.2 <- unlist(gregexpr(".nc",testForPeriod))
        testForPeriodLoc.2 <- testForPeriodLoc.2[length(testForPeriodLoc.2)]
        
        testForPeriodLoc.1 <-  as.numeric(substr(testForPeriod,testForPeriodLoc.1+1,testForPeriodLoc.1+4))
        testForPeriodLoc.2 <- as.numeric(substr(testForPeriod,testForPeriodLoc.2-6,testForPeriodLoc.2-3))
        
        if( is.na(testForPeriodLoc.1) | is.na(testForPeriodLoc.2) | (testForPeriodLoc.1 > testForPeriodLoc.2) ) { next }
        
        testForPeriod <-testForPeriodLoc.1:testForPeriodLoc.2
        
        if( ! TRUE %in% (testForPeriod %in% period)) { next }
        
        nameFile <- c(gregexpr(comb[f,"variable"],wgetFile.i)[[1]][1], ifelse(gregexpr("\\.nc",wgetFile.i)[[1]][1] < 10 , gregexpr("\\.nc",wgetFile.i)[[1]][2] , gregexpr("\\.nc",wgetFile.i)[[1]][1]) )
        nameFile <- substr(wgetFile.i, nameFile[1], nameFile[2]-1)
        
        downloadURLS <- gregexpr("http",wgetFile.i)[[1]][1]
        downloadURLS <- c(downloadURLS,gregexpr("\\.nc'",wgetFile.i)[[1]][which(gregexpr("\\.nc'",wgetFile.i)[[1]] > downloadURLS)[1]])
        downloadURLS <- substr(wgetFile.i, downloadURLS[1], downloadURLS[2]+2)
        
        if( ! is.na(nameFile) ) { passError <- TRUE }
        if( is.na(nameFile) ) { passError <- FALSE }
        
        cat(nameFile,"\n")
        attempts <- 0
        
        while( passError ) {
          
          passError <- FALSE
          attempts <- attempts + 1
          tryCatch( download.file(downloadURLS,destfile=paste0(folder[1],"/",nameFile,".nc"), method ="curl",quiet=TRUE, timeout=1000000,extra='-L -k') , error=function(e){ passError <<- TRUE } )
          if( passError ) { tryCatch( download.file(downloadURLS,destfile=paste0(folder[1],"/",nameFile,".nc"), method ="curl",quiet=TRUE, timeout=1000000) , error=function(e){ passError <<- TRUE } ) }
          if( ! passError & file.info(paste0(folder[1],"/",nameFile,".nc"))$size < 10000 ) { passError <- TRUE; file.remove(paste0(folder[1],"/",nameFile,".nc")) }
          if( attempts > 99 ) { break }
          
        }
        
      }
      
    }
    
    return(NULL)
    
  }
  
  stopCluster(cl.2) ; rm(cl.2) ; gc(reset=TRUE)
  
}

## ---------------------------

rasterFromXYZ.V2 <- function(coords,values) {
  
  tryCatch( rawDataVarTestPlot <- rasterFromXYZ(as.matrix(cbind(coords,values))) , error=function(e){ 
    
    fact <- round(rawDataResolution / res(shape)[1])
    
    if(fact > 1) { shape.i <- aggregate(shape,fact) }  
    if(fact == 1) { shape.i <- aggregate(shape,2) }  
    if(fact < 1) { shape.i <- disaggregate(shape,round(1/(rawDataResolution / res(shape)[1]))) }
    if(fact == 0) { shape.i <- shape }
    
    rawDataVarTestPlot <- rasterize( coords,shape.i,values,fun=mean,background=NA,na.rm=TRUE) 
    
    return(raster::resample(rawDataVarTestPlot, shape.i, method="bilinear",na.rm=TRUE))
    
  } )
  
}

## ---------------------------

rotateMatrix <- function(x) t(apply(x, 2, rev))

## ---------------------------

ncCreate <- function(path, 
                     xmax = NULL, xmin = NULL, xlength = NULL,
                     ymax = NULL, ymin = NULL, ylength = NULL,
                     realm = NULL,
                     experiment = NULL,
                     period = NULL,
                     save = FALSE ){
  
  # Create file
  nc <- RNetCDF::create.nc(path, format = "netcdf4")
  
  # Global attributes
  RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "Conventions", type = "NC_CHAR", value = "CF-1.5")
  RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "title", type = "NC_CHAR", value = paste0("Bio-Oracle ",variableLongName," [", realm,"]."))
  RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "experiment", type = "NC_CHAR", value = experiment)
  RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "standard_name", type = "NC_CHAR", value = variableStandardName)
  
  RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "institution", type = "NC_CHAR", value = "Bio-Oracle consortium: https://www.bio-oracle.org")
  RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "source", type = "NC_CHAR", value = "Bio-Oracle version V3.0")
  RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "history", type = "NC_CHAR", value = paste("File created:" , Sys.time() ))
  RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "comment", type = "NC_CHAR", value = "Uses attributes recommended by http://cfconventions.org")
  RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "references", type = "NC_CHAR", value = "https://www.bio-oracle.org")
  
  # Set crs
  RNetCDF::var.def.nc(nc, varname = "crs", vartype = "NC_CHAR", dimensions = NA)
  RNetCDF::att.put.nc(nc, variable = "crs", name = "grid_mapping_name", type = "NC_CHAR", value = "latitude_longitude")
  RNetCDF::att.put.nc(nc, variable = "crs", name = "long_name", type = "NC_CHAR", value = "CRS definition")
  RNetCDF::att.put.nc(nc, variable = "crs", name = "longitude_of_prime_meridian", type = "NC_DOUBLE", value = 0.)
  RNetCDF::att.put.nc(nc, variable = "crs", name = "semi_major_axis", type = "NC_DOUBLE", value = 6378137.)
  RNetCDF::att.put.nc(nc, variable = "crs", name = "inverse_flattening", type = "NC_DOUBLE", value = 298.257223563)
  RNetCDF::att.put.nc(nc, variable = "crs", name = "spatial_ref", type = "NC_CHAR", value = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]')
  RNetCDF::att.put.nc(nc, variable = "crs", name = "GeoTransform", type = "NC_CHAR", value = '-180 0.08333333333333333 0 90 0 -0.08333333333333333 ')
  
  # Set Longitude dimension
  RNetCDF::dim.def.nc(nc, dimname = "lon", dimlength = xlength)
  RNetCDF::var.def.nc(nc, varname = "lon", vartype = "NC_FLOAT", dimensions = "lon")
  RNetCDF::att.put.nc(nc, variable = "lon", name = "units", type = "NC_CHAR", value = "degrees_east")
  RNetCDF::att.put.nc(nc, variable = "lon", name = "long_name", type = "NC_CHAR", value = "Longitude")
  RNetCDF::att.put.nc(nc, variable = "lon", name = "standard_name", type = "NC_CHAR", value = "longitude")
  RNetCDF::att.put.nc(nc, variable = "lon", name = "axis", type = "NC_CHAR", value = "X")
  RNetCDF::att.put.nc(nc, variable = "lon", name = "reference_datum", type = "NC_CHAR", value = "geographical coordinates, WGS84 projection")
  RNetCDF::att.put.nc(nc, variable = "lon", name = "valid_min", type = "NC_DOUBLE", value = -180.0)
  RNetCDF::att.put.nc(nc, variable = "lon", name = "valid_max", type = "NC_DOUBLE", value = 180.0)
  
  lon_data = seq(from = xmin, to = xmax, length.out = xlength)
  RNetCDF::var.put.nc(nc, variable = "lon", data = lon_data)
  
  # Set Latitude dimension
  RNetCDF::dim.def.nc(nc, dimname = "lat", dimlength = ylength)
  RNetCDF::var.def.nc(nc, varname = "lat", vartype = "NC_FLOAT", dimensions = "lat")
  RNetCDF::att.put.nc(nc, variable = "lat", name = "units", type = "NC_CHAR", value = "degrees_north")
  RNetCDF::att.put.nc(nc, variable = "lat", name = "long_name", type = "NC_CHAR", value = "Latitude")
  RNetCDF::att.put.nc(nc, variable = "lat", name = "standard_name", type = "NC_CHAR", value = "latitude")
  RNetCDF::att.put.nc(nc, variable = "lat", name = "axis", type = "NC_CHAR", value = "Y")
  RNetCDF::att.put.nc(nc, variable = "lat", name = "reference_datum", type = "NC_CHAR", value = "geographical coordinates, WGS84 projection")
  RNetCDF::att.put.nc(nc, variable = "lat", name = "valid_min", type = "NC_DOUBLE", value = -90.0)
  RNetCDF::att.put.nc(nc, variable = "lat", name = "valid_max", type = "NC_DOUBLE", value = 90.0)
  
  lat_data = seq(from = ymin, to = ymax, length.out = ylength)
  RNetCDF::var.put.nc(nc, variable = "lat", data = lat_data)
  
  # Set time dimension
  RNetCDF::dim.def.nc(nc, dimname = "time", dimlength = length(period)) # 3 years levels
  RNetCDF::var.def.nc(nc, varname = "time", vartype = "NC_INT", dimensions = "time")
  RNetCDF::att.put.nc(nc, variable = "time", name = "units", type = "NC_CHAR", value = "days since 1970-01-01 00:00:00")
  RNetCDF::att.put.nc(nc, variable = "time", name = "long_name", type = "NC_CHAR", value = "time")
  RNetCDF::att.put.nc(nc, variable = "time", name = "calendar", type = "NC_CHAR", value = "gregorian")
  RNetCDF::att.put.nc(nc, variable = "time", name = "standard_name", type = "NC_CHAR", value = "time")
  RNetCDF::var.put.nc(nc, variable = "time", data = period)
  
  cat(paste0("# .nc file saved to disk \n"))

  return(nc)
}

## ---------------------------

ncAddVar <- function(nc, varname, varUnits, long_name, crs, compression=NA){
  RNetCDF::var.def.nc(nc, varname = varname, vartype = "NC_DOUBLE", dimensions = c("lon", "lat", "time"), deflate=compression)
  RNetCDF::att.put.nc(nc, variable = varname, name = "_FillValue", type = "NC_DOUBLE", value = -9999.9)
  RNetCDF::att.put.nc(nc, variable = varname, name = "units", type = "NC_CHAR", value = varUnits)
  RNetCDF::att.put.nc(nc, variable = varname, name = "long_name", type = "NC_CHAR", value = long_name)
  # RNetCDF::att.put.nc(nc, variable = varname, name = "standard_name", type = "NC_CHAR", value = standard_name)
  RNetCDF::att.put.nc(nc, variable = varname, name = "coordinates", type = "NC_CHAR", value = "lat lon")
  RNetCDF::att.put.nc(nc, variable = varname, name = "grid_mapping", type = "NC_CHAR", value = "crs")
  sync.nc(nc)
}

## ---------------------------

getDimNames <- function(nc.file){
  
  if( exists("dimNameLon") ) { rm(dimNameLon, pos = ".GlobalEnv"); rm(dimNameLat, pos = ".GlobalEnv"); rm(dimNameDepth, pos = ".GlobalEnv") }

  nc <- nc_open( nc.file, readunlim=FALSE )
  nc.dimension.names <- NcReadVarNames(nc)
  
  if ("height" %in% nc.dimension.names ) { dimNameDepth <<- "height" }
  if ("depth" %in% nc.dimension.names ) { dimNameDepth <<- "depth" }
  if ("lev" %in% nc.dimension.names ) { dimNameDepth <<- "lev" }
  if ("olevel" %in% nc.dimension.names ) { dimNameDepth <<- "olevel" }

  if ("longitude" %in% nc.dimension.names ) { dimNameLon <<- "longitude" }
  if ("latitude" %in% nc.dimension.names ) { dimNameLat <<- "latitude" }
  
  if ("lon" %in% nc.dimension.names ) { dimNameLon <<- "lon" }
  if ("lat" %in% nc.dimension.names ) { dimNameLat <<- "lat" }

  if ("nav_lon" %in% nc.dimension.names ) { dimNameLon <<- "nav_lon" }
  if ("nav_lat" %in% nc.dimension.names ) { dimNameLat <<- "nav_lat" }
  
  if( ! exists("dimNameDepth") ) { dimNameDepth <<- "NA" }
  
  if ( ! exists("dimNameLon") | ! exists("dimNameLat") | ! exists("dimNameDepth") ) { stop("Error :: Undefined dim names [getDimNames]") }
  
  nc_close(nc) 
  
}

## ---------------------------

getDimSizes <- function(nc.file){
  
  if( exists("dim.i") ) { rm(dim.i, pos = ".GlobalEnv"); rm(dim.j, pos = ".GlobalEnv"); rm(dim.z, pos = ".GlobalEnv"); rm(dim.t, pos = ".GlobalEnv") }
  
  raw.data <- nc_open( nc.file, readunlim=FALSE )
  nc.dimensions <- raw.data$dim
  nc.dimension.size.names <- names(raw.data$dim)
  nc_close( raw.data )
  
  if ("lon" %in% nc.dimension.size.names ) { dimName.i <- "lon" }
  if ("lat" %in% nc.dimension.size.names ) { dimName.j <- "lat" }
  
  if ("nlon" %in% nc.dimension.size.names ) { dimName.i <- "nlon" }
  if ("nlat" %in% nc.dimension.size.names ) { dimName.j <- "nlat" }
  
  if ("i" %in% nc.dimension.size.names ) { dimName.i <- "i" }
  if ("j" %in% nc.dimension.size.names ) { dimName.j <- "j" }
  
  if ("ni" %in% nc.dimension.size.names ) { dimName.i <- "ni" }
  if ("nj" %in% nc.dimension.size.names ) { dimName.j <- "nj" }
  
  if ("longitude" %in% nc.dimension.size.names ) { dimName.i <- "longitude" }
  if ("latitude" %in% nc.dimension.size.names ) { dimName.j <- "latitude" }
  
  if ("x" %in% nc.dimension.size.names ) { dimName.i <- "x" }
  if ("y" %in% nc.dimension.size.names ) { dimName.j <- "y" }
  
  if ("time" %in% nc.dimension.size.names ) { dimNameTime <- "time" } else { dimNameTime <- NULL }

  dim.i <<- nc.dimensions[dimName.i][[1]]$len
  dim.j <<- nc.dimensions[dimName.j][[1]]$len
  dim.z <<- nc.dimensions[dimNameDepth][[1]]$len
  
  if(!is.null(dimNameTime)) { dim.t <<- nc.dimensions[dimNameTime][[1]]$len } else { dim.t <- NULL }

  if( ! rawDataIs3D | (length(depthsToInterpolate) == 1 & "depthSurf" %in% depthsToInterpolate) ) { dim.z <<- 1 }
  
  if ( is.null(dim.i) | is.null(dim.j) | is.null(dim.z) | ! exists("dim.i") | ! exists("dim.j") | ! exists("dim.z") ) { stop("Error :: Undefined dim sizes [getDimSizes]") }
  
}

## ---------------------------

getDimUnits <- function(nc.file,dimName){
  
  raw.data <- nc_open( nc.file, readunlim=FALSE )
  DimUnits <- ncatt_get(raw.data,dimName,"units")$value
  nc_close( raw.data )
  return(DimUnits)
  
}

## ---------------------------

rotateMatrixCC <- function(x) { apply(t(x),2,rev) }

## ---------------------------

ncToRaster <- function(ncFile, layer, period = period.n, rasterFile){
  
  dataNC <- nc_open( ncFile , readunlim=FALSE )
  dataVar <- ncvar_get( dataNC, layer)
  lonsDataVar <- ncvar_get( dataNC, "lon") 
  latsDataVar <- ncvar_get( dataNC, "lat")
  timeDataVar <- ncvar_get( dataNC, "time")
  nc_close( dataNC )
  
  if( is.na( dim(dataVar)[3]) ) { rasterVar <- raster(rotateMatrixCC(dataVar)) }
  if( ! is.na( dim(dataVar)[3]) ) { rasterVar <- raster(rotateMatrixCC(dataVar[,,period])) }
  
  extent(rasterVar) <- c(min(lonsDataVar), max(lonsDataVar), min(latsDataVar), max(latsDataVar))
  crs(rasterVar) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  return(rasterVar)
  
}

## ------------------------------
## ------------------------------

cat("\n")
cat("# Reading functions and environment variables","\n")
