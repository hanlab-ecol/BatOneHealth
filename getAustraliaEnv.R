#'This function operates to grab a single raster layer or multiple rasters for
#' a given variable
#' @param path This is the path to a repository with all of the environmental data used in this project. By default, the function will grab the current working directory.
#' @param variable This is the variable that you wish to query. There are a handful of options, but if you stray from one of these, the function will provide an error (because these data don't exist!)
#' @param lag An argument specifying what type of lag you wish to query. By default 
#' it is set to NA which will specify no lag, but options are typically 2, 3, 6, 9, 12,
#' 15, 18, 21, or 24 except in the case of NDVI which has options for 1, 2, 3, and 4,
#' months.
#' @param year Specifies which year of data you wish to access. For most variables, 
#' the available options are 1996 - 2021, except for land cover, which is only available
#' for the years 2002 - 2014.
#' @param month Specifies the month to grab the data from. It is set as NULL by default
#' which means that this function will grab all months of the year specified. 
#' @return This function will return an object of terra::SpatRaster-class ranging from 
#' 1 to 12 layers representing the variable chosen at the time points (year and 
#' possibly month) chosen. 
#' @example getAustraliaEnv(path = "~/Desktop/AustralianClimateData", variable = "TemperatureMax", lag = 6, year = 2019, month = 8)
getAustraliaEnv <- function(path = getwd(), variable = c("TemperatureMax", "TemperatureMin", "TemperatureDiff", "Precipitation", "PrecipitationAnomaly", "NDVI", "VaporPressure9AM", "SolarExposure", "SoilMoistureRootZone", "Evapotranspiration", "LandCover"), lag = NA, year, month = NULL) {
  if(!variable %in% c("TemperatureMax", "TemperatureMin", "TemperatureDiff", "Precipitation", "PrecipitationAnomaly", "NDVI", "VaporPressure9AM", "SolarExposure", "SoilMoistureRootZone", "Evapotranspiration", "LandCover")) {
    stop("Variable name, `", variable, "`, is not available to attach. Please check spelling and list of available variables.")
  }
  PATH <- paste0(path, variable, "Lag")
  if(variable %in% c("SoilMoistureRootZone", "Evapotranspiration", "LandCover")) {
    PATH <- paste0(path, variable)
  }
  if(!is.na(lag)) {
    PATH <- paste0(PATH, "/", variable, "Lag", lag, "Months")
  }
  if(is.null(month)) {
    terra::rast(list.files(path = PATH, pattern = as.character(year), full.names = T))
  } else {
    terra::rast(list.files(path = PATH, pattern = as.character(year), full.names = T))[[month]]
  }
}