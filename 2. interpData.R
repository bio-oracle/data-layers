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
cat("# Interpolation //",depthInterpolation,"\n")
time.i <- Sys.time()

## ------------------------
## ------------------------

seqListing <- raw.data.depth
seqChunks <- data.frame(from = c(0,seqListing[-length(seqListing)]), to = c(seqListing[-length(seqListing)],11000) )

if( realm == "Surface") { seqChunks <- seqChunks[1,] }

## ------------------------

rawDataCellsDT <- as.data.table(interpSurfaceMetric)
interpCellsDT <- interpCells[,.(Lon,Lat,get(depthInterpolation))]
colnames(interpCellsDT) <- c("Lon","Lat","Depth")

## ------------------------

for( d in 1:nrow(seqChunks) ) {
        
    interpCellsDT.d.index <- interpCellsDT[ Depth >= seqChunks[d,1] & Depth < seqChunks[d,2] , which = TRUE  ]
    interpCellsDT.d <- interpCellsDT[ interpCellsDT.d.index , .(Lon,Lat,Depth)  ]
    interpCellsDT.d[, Lon := Lon*1/0.00001]
    interpCellsDT.d[, Lat := Lat*1/0.00001]
    
    if( nrow(interpCellsDT.d ) == 0 ) { next }
    
    ## ----------------------
    
    rawDataCellsDT.d <- rawDataCellsDT[ Depth >= seqChunks[d,1] - 0.001 & Depth <= seqChunks[d,2] + 0.001 ,]
    
    depthChunkCorrection <- 0
    while( nrow(rawDataCellsDT.d) == 0 ) {
      depthChunkCorrection <- depthChunkCorrection + 1
      rawDataCellsDT.d <- rawDataCellsDT[ Depth >= seqChunks[max(1,seqChunks[d,1] - depthChunkCorrection),1] - 0.001 & Depth < seqChunks[max(1,seqChunks[d,2] - depthChunkCorrection),2] + 0.001 ,]
    }
    
    rawDataCellsDT.d[, Lon := Lon*1/0.00001]
    rawDataCellsDT.d[, Lat := Lat*1/0.00001]
    setkey(rawDataCellsDT, Lon,Lat)
    
    # ----------------------------------
    
    seqListing <- seq(1,nrow(interpCellsDT.d),by=5000)
    if( length(seqListing) < nCores) { seqListing <- round(seq(1,nrow(interpCellsDT.d),length.out=32)) }
    
    seqListing[length(seqListing)] <- nrow(interpCellsDT.d)
    parallelChunks <- data.frame(from = seqListing[-length(seqListing)], to = c(seqListing[-c(1,length(seqListing))] - 1 , nrow(interpCellsDT.d) ) )
    
    Cluster <- makeCluster( nCores )
    registerDoParallel(Cluster)
    
    interpProcessKN <- foreach(parallelChunk=1:nrow(parallelChunks), .verbose=FALSE, .packages=c("ncdf4","data.table","bigmemory","FNN")) %dopar% {
            
          cellToInterpChunk <- parallelChunks[parallelChunk,1]:parallelChunks[parallelChunk,2]
          cellToInterpChunk <- unique(cellToInterpChunk[cellToInterpChunk != 0])
          
          idw.nearest.r <- get.knnx( rawDataCellsDT.d[ , .(Lon,Lat) ] , interpCellsDT.d[ cellToInterpChunk, .(Lon,Lat)], k= ifelse( realm == "Surface" , (interpNMax / 2) - 2 , interpNMax - 4 ) , algorithm="kd_tree" )
          idw.nearest.i <- idw.nearest.r$nn.index
          idw.nearest.d <- idw.nearest.r$nn.dist
    
          for( pred in outcomePeriodPredictors) {   
            assign(paste0("rawDataVar",pred,"Df"),data.frame(matrix(NA,nrow=length(cellToInterpChunk),ncol=dim.t)))
          }
          
          raw.data <- nc_open( paste0(tempFolder,"/","rawDataDepth_",d,".nc"), readunlim=FALSE )
    
          for(i in 1:length(cellToInterpChunk)) {
            
              idw.nearest <- rawDataCellsDT.d[ idw.nearest.i[i,], .(Var1 , Var2, Var3 ) ]
              idw.nearest[ , Var3 := Var3 - min(Var3) + 1 ]
              
              if( realm != "Surface" ) {
                
                verticaldistance <- abs(as.numeric(interpCellsDT.d[ cellToInterpChunk[i], .(Depth)]) - seqChunks[d,])
                idw.nearest.d[i,1:ncol(idw.nearest.d)] <- as.numeric(idw.nearest.d[i,1:ncol(idw.nearest.d)] + verticaldistance[idw.nearest$Var3])
                
              }
              
              idw.nearest.d[idw.nearest.d == 0] <- 0.0000000001

              rawDataVarMin <- apply(idw.nearest,1,function(x) { ncvar_get( raw.data, "raw.dataMin", start=c(x[1],x[2],x[3],1), count=c(1,1,1,dim.t)) })
              rawDataVarMinDf[i,] <- sapply(1:dim.t, function(f) { (sum( rawDataVarMin[f,] / idw.nearest.d[i,]^idwPower , na.rm=TRUE)) / (sum( 1 / idw.nearest.d[i,]^idwPower , na.rm=TRUE)) } )  
              
              rawDataVarMax <- apply(idw.nearest,1,function(x) { ncvar_get( raw.data, "raw.dataMax", start=c(x[1],x[2],x[3],1), count=c(1,1,1,dim.t)) })
              rawDataVarMaxDf[i,] <- sapply(1:dim.t, function(f) { (sum( rawDataVarMax[f,] / idw.nearest.d[i,]^idwPower , na.rm=TRUE)) / (sum( 1 / idw.nearest.d[i,]^idwPower , na.rm=TRUE)) } )  
              
              rawDataVarMean <- apply(idw.nearest,1,function(x) { ncvar_get( raw.data, "raw.dataMean", start=c(x[1],x[2],x[3],1), count=c(1,1,1,dim.t)) })
              rawDataVarMeanDf[i,] <- sapply(1:dim.t, function(f) { (sum( rawDataVarMean[f,] / idw.nearest.d[i,]^idwPower , na.rm=TRUE)) / (sum( 1 / idw.nearest.d[i,]^idwPower , na.rm=TRUE)) } )  
    
              if( "ltMax" %in% outcomePeriodPredictors & "ltMin" %in% outcomePeriodPredictors ) {
                
                rawDataVarltMax <- apply(idw.nearest,1,function(x) { ncvar_get( raw.data, "raw.dataLtMax", start=c(x[1],x[2],x[3],1), count=c(1,1,1,dim.t)) })
                rawDataVarltMaxDf[i,] <- sapply(1:dim.t, function(f) { (sum( rawDataVarltMax[f,] / idw.nearest.d[i,]^idwPower , na.rm=TRUE)) / (sum( 1 / idw.nearest.d[i,]^idwPower , na.rm=TRUE)) } )  
                
                rawDataVarltMin <- apply(idw.nearest,1,function(x) { ncvar_get( raw.data, "raw.dataLtMin", start=c(x[1],x[2],x[3],1), count=c(1,1,1,dim.t)) })
                rawDataVarltMinDf[i,] <- sapply(1:dim.t, function(f) { (sum( rawDataVarltMin[f,] / idw.nearest.d[i,]^idwPower , na.rm=TRUE)) / (sum( 1 / idw.nearest.d[i,]^idwPower , na.rm=TRUE)) } )  
            
                rawDataVarRangeDf[i,] <- rawDataVarMaxDf[i,] - rawDataVarMinDf[i,]
                
              }
          }
          
          nc_close( raw.data )
          
          ## --------------------
          
          for( pred in outcomePeriodPredictors) {   
            
            climatologyDescFile <- paste0(processTempFolder,"/climatology",pred,depthInterpolation,".desc")
            assign("climatologyDesc",dget( climatologyDescFile) )
            climatologyDescBM <- attach.big.matrix( climatologyDesc )

            for(t.i in 1:dim.t) {  
              
              valesToDump <- get(paste0("rawDataVar",pred,"Df"))[,t.i]
              
              if( variable == "KDPAR_mean" | variable == "PAR_mean" ) {
                valesToDump[valesToDump == 0] <- NA
              }
              
              climatologyDescBM[interpCellsDT.d.index[cellToInterpChunk],t.i] <- valesToDump
    
          } }
          
          ## ------------
          
          for( pred in outcomePeriodPredictors) {   
            rm(list = paste0("rawDataVar",pred,"Df"))
            rm(list = paste0("rawDataVar",pred))
          }
          rm(climatologyDescBM)
          
          gc(reset=TRUE)
          return(NULL)
          
    }
    
    stopCluster(Cluster) ; rm(Cluster)
    
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

shape <- raster(paste0(bathymetryFolder,"/",bathymetryFiles[1]))
shape <- crop(shape,extent(outcomeExtent))
shape[!is.na(shape)] <- 1 

## ------------------------------
## ------------------------------

climatologyMean <- dget( paste0(processTempFolder,"/","climatology","Mean",depthInterpolation,".desc"))
climatologyMean <- attach.big.matrix(climatologyMean)
climatologyCells <- dget( paste0(processTempFolder,"/","cellsStructure.desc"))
climatologyCells <- attach.big.matrix(climatologyCells)
climatologyCells <- bigmemory::as.matrix(climatologyCells)
colnames(climatologyCells) <- c("Lon","Lat","Row","Col","Cell","depthSurf","depthMin","depthMean","depthMax")

shape[climatologyCells[,5]] <- climatologyMean[,ncol(climatologyMean)]
plot(shape)

png(filename = paste0(processTempFolder,"/interpDataTestImage",depthInterpolation,".png"), width = 1920, height = 1080)
plot(shape)
dev.off()
writeRaster(shape,filename=paste0(processTempFolder,"/interpDataTestImage",depthInterpolation,".tif"),format="GTiff",overwrite=T)

closeAllConnections()
gc(reset=TRUE, full = TRUE)

## ------------------ 

cat("# End process // ET:", difftime(Sys.time(), time.i, units='hours') ,"hours\n")

## --------------------------------------------------------
## --------------------------------------------------------