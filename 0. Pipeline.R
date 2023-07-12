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

## Git
## library(credentials)
## set_github_pat()

## ------------------------------------------------------------------ ##
## ------------------------------------------------------------------ ##

closeAllConnections()
rm(list=ls())
gc(reset=TRUE)

useArmadillo <- TRUE

source("Dependencies/mainFunctions.R")
source("0. mainConfiguration.R")

nCores <- 10

## ------------------------------ ##
## Main configuration

exportFolder
mainResolution
bathymetryFiles[2]

outcomePeriodPredictors <- c("Max","ltMax","Mean","ltMin","Min","Range") # c("Max","ltMax","Mean","ltMin","Min","Range") // c("Max","Mean","Min")

outcomePeriodY <- 2000:2020 # 2000:2020 // 2020:2100
outcomePeriodM <- 1:12
outcomePeriodMName <- NULL # "seasonName"
exportClimatologyRasters <- FALSE
baselinePeriod <- 2010:2020

## ------------------------------ ##

# listExperiments$name
# listModels$name
# listModelsFinal$name
# listVariables$name
# listDataSources$name

pipeLineLoop <- expand.grid( variable = c("KDPAR_mean") , # PAR_mean KDPAR_mean
                             dataSource = c("GlobColour") , # "myOcean" // "ecmwf" // "cmip6" // "GlobColour"
                             periodType = "decade", # "year" // "decade"
                             experiment =  listExperiments[c(4),1] , #  //"Baseline" listExperiments[c(10:15),1]
                             model = NA , # NA listModelsFinal$name
                             realm = c("Surface"), # Benthic Surface Both
                             stringsAsFactors = FALSE)

pipeLineLoop <- pipeLineLoop[ sort(pipeLineLoop$variable,index.return=T)$ix, ]
pipeLineLoop

## ------------------------------ ##
## If CMIP data

availableRawData <- ListingAvailableRawData(folder = paste0(rawDataFolder,"/CMIP6"))
pipeLineLoop <- pipeLineLoop[which(sapply(which(grepl("cmip",pipeLineLoop$dataSource)),function(x) { ifelse(length(which(availableRawData$variable == pipeLineLoop[x,"variable"] & availableRawData$experiment == pipeLineLoop[x,"experiment"] & availableRawData$model == pipeLineLoop[x,"model"])) > 0 , TRUE, FALSE) } )),]
pipeLineLoop

## --------

unique(availableRawData$dimNameDepth)
unique(availableRawData$dimNameGeo)

selectModels <- aggregate(pipeLineLoop$experiment,list(pipeLineLoop$model),length); selectModels
selectModels <- selectModels[selectModels$x >= 3,1]; selectModels
pipeLineLoop <- pipeLineLoop[pipeLineLoop$model %in% selectModels,]
pipeLineLoop

pipeLineLoop <- pipeLineLoop[!pipeLineLoop$experiment %in% c("ssp119","ssp585"),]
pipeLineLoop

## ------------------------------ ##

# unlink(tempFolder, recursive = TRUE)
# rm(rawDataReady)
# closeAllConnections()
# gc(reset=TRUE)

overwrite <- TRUE

for( pipe.i in 1:nrow(pipeLineLoop)) {

  ## ---------------- ##
  
  preparingClimatologies <- TRUE
  
  ## ---------------- ##
  
  variable <- pipeLineLoop[pipe.i,"variable"]
  dataSource <- pipeLineLoop[pipe.i,"dataSource"]
  model <- pipeLineLoop[pipe.i,"model"]
  experiment <- pipeLineLoop[pipe.i,"experiment"]
  realm.i <- pipeLineLoop[pipe.i,"realm"]
  outcomePeriodType <- pipeLineLoop[pipe.i,"periodType"]
  
  ## ---------------- ##
  
  cat(paste0( pipe.i , " out of " , nrow(pipeLineLoop), " :: ", variable, " | " , dataSource, " | " , model, " | " , experiment, " | " , realm.i, " | " , outcomePeriodType ), file = "logFile.txt", sep = "\n", append = TRUE)
  
  ## ---------------- ##
  
  if( pipe.i == 1 | ! exists("startTime") ) { startTime <- Sys.time() }

  cat('\014')
  cat('\n')
  cat("# ---------------")
  cat('\n')
  cat('# Process:',pipe.i,"out of",nrow(pipeLineLoop),'\n')
  cat('# Time taken:',as.numeric(Sys.time() - startTime),'\n')
  cat('# Variable:',variable,"\n")
  cat('# Experiment:',experiment,"\n")
  cat('# Model:' ,ifelse(!is.na(model),model,"NA"))
  cat("\n")
  cat("# ---------------")
  cat("\n")
  cat("#", baselinePeriod[length(baselinePeriod)],"->",outcomePeriodY[length(outcomePeriodY)])
  cat("\n")
  cat("# ---------------")
  cat("\n")
  cat("\n")
  
  ## ---------------- ##
  
  depthsToInterpolate <- character(0)
  if( realm.i == "Surface" | realm.i == "Both" ) { depthsToInterpolate <- "depthSurf" } 
  if( realm.i == "Benthic" | realm.i == "Both" ) { depthsToInterpolate <- c(depthsToInterpolate,"depthMin","depthMean","depthMax") } 
  
  ## ---------------- ##
  
  processedData <- list.files(paste0(exportFolders,"/mainClimatologies"), full.names = TRUE)
  processedData <- processedData[grepl(variable,processedData)]
  processedData <- processedData[grepl(experiment,processedData)]
  processedData <- list.files(processedData, recursive=TRUE, full.names = TRUE, pattern="\\.nc")
  processedData <- processedData[grepl(model,processedData)]
  processedData <- sapply(depthsToInterpolate, function(x) TRUE %in% grepl(firstToupper(x),processedData) )
  
  if( length(unlist(processedData)) != 0 ) {  if( sum(processedData) == length(depthsToInterpolate) & ! overwrite ) { next } }
  
  ## ---------------- ##
  
  for( depthInterpolation in depthsToInterpolate ) {
    
    if( depthInterpolation == "depthSurf" ) { realm <- "Surface" } 
    if( depthInterpolation %in% c("depthMin","depthMean","depthMax") ) { realm <- "Benthic" } 
    
    source("1. defineEnvironment.R")
    
    if( ! exists("rawDataReady") ) { source("1. prepareRawData.R") }
    
    if( depthInterpolation %in% c("depthMin","depthMean","depthMax") & ! rawDataIs3D ) { next }
     
    source("1. prepareClimatology.R")
    source("2. interpData.R")
    
    ## ---------------- ##

    cat("# Moving files to final path \n")
    
    exportFolderClimatologies <- paste0(exportFolder,"/mainClimatologies/",climatologyName,ifelse( experiment != "Baseline" & (model != "NA" | ! is.na(model) ) , paste0("/Models/",model) , "" ),"/")
    if( depthInterpolation == depthsToInterpolate[1] & dir.exists( exportFolderClimatologies ) ) { unlink(exportFolderClimatologies, recursive = TRUE) }
    
    attempt <- 0
    while( length(list.files(exportFolderClimatologies, pattern=paste0(depthInterpolation,".bin"))) == 0 ) {
      attempt <- attempt + 1
      if( attempt == 99 ) { stop("Error :: 111") }
      if( ! dir.exists( exportFolderClimatologies ) ) { dir.create( exportFolderClimatologies , recursive = TRUE) }
      for( x in list.files(processTempFolder, full.names = T) ) {
        file.copy(x,paste0(exportFolderClimatologies,"/",gsub(processTempFolder,"",x)),overwrite=TRUE)
    } }

    ## ---------------- ##
    
    source("3. exportClimatologies.R")
    
  }
  
  ## ---------------- ##
  
  unlink(tempFolder, recursive = TRUE)
  rm(rawDataReady)
  
  # delete temp BM files
  
  file.remove(list.files(exportFolderClimatologies, full.names = T, pattern=".bin"))
  file.remove(list.files(exportFolderClimatologies, full.names = T, pattern=".desc"))
    
  ## ---------------- ##
  
}

## ------------------------------ ##
## ------------------------------ ##

## Ensemble climatologies :: CMIP data

pipeLineLoop <- unique(pipeLineLoop[,c("variable","dataSource","experiment","realm","periodType")])

for(pipe.i in 1:nrow(pipeLineLoop)) {
  
  ## ---------------- ##
  
  preparingClimatologies <- FALSE
  
  ## ---------------- ##
  
  variable <- pipeLineLoop[pipe.i,"variable"]
  experiment <- pipeLineLoop[pipe.i,"experiment"]
  realm.i <- pipeLineLoop[pipe.i,"realm"]
  outcomePeriodType <- pipeLineLoop[pipe.i,"periodType"]
  
  availableClimatologies <- list.files( paste0(exportFolders,"/mainClimatologies/"))
  availableClimatologies <- availableClimatologies[grepl(paste0(variable,"_"),availableClimatologies)]
  availableClimatologies <- availableClimatologies[grepl(experiment,availableClimatologies)]
  availableClimatologies <- availableClimatologies[grepl(firstToupper(outcomePeriodType),availableClimatologies)]
  
  if(length(availableClimatologies) > 1) { stop("Review :: 703") }
  
  availableClimatologiesModels <- list.dirs( paste0(exportFolders,"/mainClimatologies/",availableClimatologies,"/Models/"), full.names = TRUE )
  availableClimatologiesModels <- availableClimatologiesModels[which(apply(sapply(listModels$name,function(x) grepl(x,availableClimatologiesModels) ),1,sum) != 0)]
  
  if(length(availableClimatologiesModels) == 0) { stop("Review :: 903") }
  
  variable.m <- variable
  
  mainFolderBaselineClimatology <- list.files(paste0(exportFolders,"/mainClimatologies/"))
  mainFolderBaselineClimatology <- mainFolderBaselineClimatology[grepl("Baseline",mainFolderBaselineClimatology)]
  mainFolderBaselineClimatology <- mainFolderBaselineClimatology[grepl(paste0(variable.m,"_"),mainFolderBaselineClimatology)]
  mainFolderBaselineClimatology <- mainFolderBaselineClimatology[grepl(firstToupper(outcomePeriodType),mainFolderBaselineClimatology)]
  mainFolderBaselineClimatology <- paste0(exportFolders,"/mainClimatologies/",mainFolderBaselineClimatology)
  mainFolderBaselineClimatology <- mainFolderBaselineClimatology[ unlist(sapply(mainFolderBaselineClimatology, function(x) {  length(list.files(x)) >= 1 })) ]
  
  if(length(mainFolderBaselineClimatology) == 0) { stop("Review :: 904") }
  if(length(mainFolderBaselineClimatology) > 1) { stop("Review :: 924") }
  
  ## ---------------- ##
  
  climatologyis3D <- (TRUE %in% grepl("DepthMean.nc",list.files(availableClimatologiesModels))) & (TRUE %in% grepl("DepthMean.nc",list.files(mainFolderBaselineClimatology)))
  
  ## ---------------- ##
  
  depthsToInterpolate <- character(0)
  if( realm.i == "Surface" | realm.i == "Both" ) { depthsToInterpolate <- "depthSurf" } 
  if( realm.i == "Benthic" | realm.i == "Both" ) { depthsToInterpolate <- c(depthsToInterpolate,"depthMin","depthMean","depthMax") } 
  
  for( depthInterpolation in depthsToInterpolate ) {
    
    if( depthInterpolation == "depthSurf" ) { realm <- "Surface" } 
    if( depthInterpolation %in% c("depthMin","depthMean","depthMax") ) { realm <- "Benthic" } 
    if( depthInterpolation %in% c("depthMin","depthMean","depthMax") & ! climatologyis3D ) { next } 
    
    exportFolder <- exportFolders[which.max(apply(sapply(exportFolders, function(x) { grepl(x,availableClimatologiesModels) } ),2,sum))]
    
    source("1. defineEnvironment.R")
    source("5. ensembleClimatologies.R")
    
    exportFolder <-exportFolders[1]
    
  }
  
  # Delete model files, leave metadata and images
  
  for( x in availableClimatologiesModels ) { file.remove(list.files(x, pattern=".nc", full.names=T)) }
  
}

## ------------------------------------------------------------------ ##
## ------------------------------------------------------------------ ##
## END OF CODE