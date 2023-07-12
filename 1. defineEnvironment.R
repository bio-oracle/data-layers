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
cat("# Defining environment","\n")
time.i <- Sys.time()

## ------------------------
## ------------------------

if( experiment == "Baseline" & sum(2000:2020 %in% outcomePeriodY) < 5 ) { stop("Error :: 111")}
if( experiment != "Baseline" & sum(2000:2020 %in% outcomePeriodY) > 10 ) { stop("Error :: 112")}

## ------------------------

if(experiment == "lgm") { stop("Generate new project. Pay attention to differneces in depth ranges (-120) for benthic interpolations") }

## ------------------------

if(outcomePeriodType == "year" & "ltMax" %in% outcomePeriodPredictors) { stop("Error :: 310")}

## ------------------------

if( ! TRUE %in% grepl(gsub("\\.","",as.character(mainResolution)),bathymetryFiles[1]) ) { stop("Error :: 919")}

## ------------------------

variableLongName <<- as.character(listVariables[listVariables$name == variable,"longName"])
variableStandardName <<- as.character(listVariables[listVariables$name == variable,"standardName"])
variableMinPossibleValue <<- as.numeric(listVariables[listVariables$name == variable,"minPossibleValue"])
variableMaxPossibleValue <<- as.numeric(listVariables[listVariables$name == variable,"maxPossibleValue"])

## ------------------------
## ------------------------

if( preparingClimatologies ) {
  
  processTempFolder <- paste0(tempFolder,"/mainClimatologies/",variable,"_",experiment,ifelse(is.null(outcomePeriodMName),"",paste0(outcomePeriodMName,"_")),"_",firstToupper(outcomePeriodType),"_",paste0(min(outcomePeriodY),"_",max(outcomePeriodY)),ifelse(!is.na(model),paste0("/",model,"/"),"/"))
  if(! dir.exists( processTempFolder )) { dir.create( processTempFolder ,recursive = TRUE ) }
  if(! dir.exists( exportFolder )) { dir.create( exportFolder ,recursive = TRUE ) }
  
  if( dataSource == "myOcean" ) { 
    
    variable.i <<- variable
    mainRawDataFolder <<- paste0( rawDataFolder,"/Copernicus") 
    if( variable %in% c("thetao","so","uo","vo","mlotst","zos","siconc","sithick") ) { mainRawDataFolder <<- paste0(mainRawDataFolder,"/GLOBAL_REANALYSIS_PHY_001_030/") }
    if( variable %in% c("chl","no3","po4","si","nppv","o2","fe","phyc","ph") ) { mainRawDataFolder <<- paste0(mainRawDataFolder,"/GLOBAL_REANALYSIS_BIO_001_029/") }
    if( variable %in% c("PP") ) { mainRawDataFolder <<- paste0(mainRawDataFolder,"/OCEANCOLOUR_GLO_CHL_L4_REP_OBSERVATIONS_009_082/") }
    
  }
  
  if( dataSource == "cmip6" ) { 
    
    mainRawDataFolder <<- paste0(rawDataFolder,"/CMIP6")
    variable.i <<- variable
    
  }
  
  if(dataSource == "cmip5" ) {  stop("Review :: 000") }
  
  if( dataSource == "GlobColour" ) { 
    
    useArmadillo <<- FALSE
    cat("# useArmadillo: set to FALSE\n")
    
    variable.i <<- variable
    mainRawDataFolder <<- paste0(rawDataFolder,"/GlobColor")
    
  }
  
  if( dataSource == "ecmwf" ) {
    
    variable.i <<- variable
    
    if( variable == "u10" | variable == "v10" ) { 
      
      mainRawDataFolder <<- paste0(rawDataFolder,"/ECMWF/Wind Speed/")
      
      }
    
    if( variable == "clt" ) {
      
      variable.i <<- "lcc"
      mainRawDataFolder <<- paste0(rawDataFolder,"/ECMWF/Cloud Cover/") 
      
    }
    
    if( variable == "tas" ) {
      
      variable.i <<- "t2m"
      mainRawDataFolder <<- paste0(rawDataFolder,"/ECMWF/Air Temperature/") 
      
    }
    
    if( variable != "u10" & variable != "v10" & variable != "clt" & variable != "tas" ) { stop("Review :: 009") }
    
  }
  
}

if( ! preparingClimatologies ) {
  
  dataSource <- "Climatologies"
  
}

## ------------------
## ------------------ 

cat("# End process // ET:", difftime(Sys.time(), time.i, units='hours') ,"hours\n")

## --------------------------------------------------------
## --------------------------------------------------------