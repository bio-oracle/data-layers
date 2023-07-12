
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
    RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "title", type = "NC_CHAR", value = paste0("Bio-Oracle ",variableLongName))
    RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "experiment", type = "NC_CHAR", value = experiment)
    
    RNetCDF::att.put.nc(nc, variable = "NC_GLOBAL", name = "institution", type = "NC_CHAR", value = "Bio-Oracle Consortium: https://www.bio-oracle.org/team.php")
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
    RNetCDF::att.put.nc(nc, variable = "time", name = "units", type = "NC_CHAR", value = "years since 0000-01-01 00:00:00")
    RNetCDF::att.put.nc(nc, variable = "time", name = "long_name", type = "NC_CHAR", value = "time")
    RNetCDF::att.put.nc(nc, variable = "time", name = "standard_name", type = "NC_CHAR", value = "time")
    RNetCDF::var.put.nc(nc, variable = "time", data = period)
    
    # Set realm dimension
    RNetCDF::dim.def.nc(nc, dimname = "realm", dimlength = length(realm))
    RNetCDF::var.def.nc(nc, varname = "realm", vartype = "NC_STRING", dimensions = "realm")
    RNetCDF::att.put.nc(nc, variable = "realm", name = "units", type = "NC_CHAR", value = "dmless")
    RNetCDF::att.put.nc(nc, variable = "realm", name = "long_name", type = "NC_CHAR", value = "Realm: depth level")
    RNetCDF::att.put.nc(nc, variable = "realm", name = "standard_name", type = "NC_CHAR", value = "realm")
    RNetCDF::var.put.nc(nc, variable = "realm", data = realm)
    
    # Save to disk if requested
    if(save){
      RNetCDF::close.nc(nc)
      message(paste0("File saved to disk at "), Sys.time(), " in: ", normalizePath(path))
      message("Open with RNetCDF::open.nc(path)")
      return()
    }
    
    return(nc)
}

ncAddVar <- function(nc, varname, long_name, standard_name, crs){
  RNetCDF::var.def.nc(nc, varname = varname, vartype = "NC_DOUBLE", dimensions = c("realm", "lon", "lat", "time"))
  RNetCDF::att.put.nc(nc, variable = varname, name = "_FillValue", type = "NC_DOUBLE", value = -9999.9)
  RNetCDF::att.put.nc(nc, variable = varname, name = "units", type = "NC_CHAR", value = "deg_C")
  RNetCDF::att.put.nc(nc, variable = varname, name = "long_name", type = "NC_CHAR", value = long_name)
  RNetCDF::att.put.nc(nc, variable = varname, name = "standard_name", type = "NC_CHAR", value = standard_name)
  RNetCDF::att.put.nc(nc, variable = varname, name = "coordinates", type = "NC_CHAR", value = "lat lon")
  RNetCDF::att.put.nc(nc, variable = varname, name = "grid_mapping", type = "NC_CHAR", value = "crs")
  sync.nc(nc)
}


