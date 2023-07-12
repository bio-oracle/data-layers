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

arma3DMean <- function(data,dim) {
  
  if( dim == 3 ) { dim.i <- 1:2 }
  result <- apply(data,dim.i,mean, na.rm=T)
  
  result[result == -Inf] <- NA
  result[result == Inf] <- NA
  return(result)
  
}

arma3DMax <- function(data,dim) {
  
  if( dim == 3 ) { dim.i <- 1:2 }
  result <- apply(data,dim.i,max, na.rm=T)
  
  result[result == -Inf] <- NA
  result[result == Inf] <- NA
  return(result)
  
}

arma3DMin <- function(data,dim) {
  
  if( dim == 3 ) { dim.i <- 1:2 }
  result <- apply(data,dim.i,min, na.rm=T)
  
  result[result == -Inf] <- NA
  result[result == Inf] <- NA
  return(result)
  
}

if(exists("useArmadillo")) {
  
  if( useArmadillo ) {
    
    cppFunction('arma::mat arma3DMean(arma::cube x, int dim) {
  arma::mat result = mean( x, dim - 1 );
  return result;
  }', depends = "RcppArmadillo")
    
    cppFunction('arma::mat arma3DMin(arma::cube x, int dim) {
  arma::mat result = min( x, dim - 1 );
  return result;
  }', depends = "RcppArmadillo")
    
    cppFunction('arma::mat arma3DMax(arma::cube x, int dim) {
  arma::mat result = max( x, dim - 1 );
  return result;
  }', depends = "RcppArmadillo")
    
  }
  
  if( ! useArmadillo ) { cat("\n","# Armadillo not loaded [!]","\n") } 
  
}
