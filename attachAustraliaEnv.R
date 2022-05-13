#'A function to attach Australian environmental data to a given set of points, years, and months
#'@param points These are the points that will have data attached to them. This argument accepts either a data.frame with two columns (representing longitude and latitude) with a coordinate reference system of GDA 2020 (EPSG:8059) or an sf object with a coordinate reference system of GDA 2020.
#'@param path This is the path to the directory that includes all of the environmental data for Australia. By default, this grabs the current working directory.
#'@param variables This is a bit of a placeholder that now only takes the option to grab all variables, but may eventually be removed or used to only grab certain variables.
#'@param year A vector of years that must match the number of points in the same order as the year information associated with each point.
#'@param month A vector of months that similar to the year argument must be given in the same order as the points.
#'@return The result of this function will be a data.frame object with the environmental information attached for each point for its given year and month. 
#'@example 
#'attachAustraliaEnv(PNTS, path = "~/Desktop/AustralianClimateData/", variables = "all", year = c(2019, 2019, 2019, 2017, 2016), month = c(5, 8, 8, 7, 2))
attachAustraliaEnv <- function(points, path = getwd(), variables = c("all", ...), year, month) {
  if(sum(class(points) %in% c("sf", "data.frame")) == 0) stop("Points argument does not seem to be of class data.frame or an sf object.")
  if(nrow(points) != length(year)) stop("Vector of years does not match with the number of points given. Check the length of years and match to the number of points.")
  if(nrow(points) != length(month)) stop("Vector of months does not match with the number of points give. Check the length of years and match to the number of points.")
  if(sum(class(points) %in% "data.frame") == 1 & sum(class(points) %in% "sf") == 0) {
    if(max(points) <= 180) message("The maximum coordinate value is less than or equal to 180, which means that these points may not be in the Geodetic Datum of Australia 2020. As a result, the environmental data may not match with these points and will result in a data.frame of mostly NAs. Please transform these data and try again.")
  }
  if(sum(class(points) %in% "sf") == 1) {
    if(st_crs(points)$input != "EPSG:8059") stop("The coordinate reference system for this object is not equivalent to the environmental data. Please transform these points using `sf::st_transform` to EPSG:8059 and then rerun.")
    points <- data.frame(st_coordinates(points))
  }
  DIST <- cbind.data.frame(year = year,
                           month = month) %>%
    distinct(year, month)
  LAGS <- list(c(NA, 2, 3, 6, 9, 12, 15, 18, 21, 24),
               c(NA, 2, 3, 6, 9, 12, 15, 18, 21, 24),
               c(NA, 2, 3, 6, 9, 12, 15, 18, 21, 24),
               c(NA, 2, 3, 6, 9, 12, 15, 18, 21, 24),
               c(NA, 2, 3, 6, 9, 12, 15, 18, 21, 24),
               c(NA, 1, 2, 3, 4),
               c(NA),
               c(NA, 2, 3, 6, 9, 12, 15, 18, 21, 24),
               c(NA),
               c(NA),
               c(NA))
  VARS <- c("TemperatureMax", "TemperatureMin", "TemperatureDiff", "Precipitation", "PrecipitationAnomaly", "NDVI", "VaporPressure9AM", "SolarExposure", "SoilMoistureRootZone", "Evapotranspiration", "LandCover")
  if(variables == "all") {
    for(m in 1:nrow(DIST)) {
      nums <- which(year == DIST[m, 1] & month == DIST[m, 2])
      pnt <- points[nums, ]
      VALS <- lapply(1:length(VARS), function(i) {
        do.call(c, lapply(LAGS[[i]], function(k) {
          if(VARS[i] == "LandCover") {
            YR <- DIST[m, 1]
            if(DIST[m, 1] > 2014) YR <- 2014
            if(DIST[m, 1] < 2002) YR <- 2002
            getAustraliaEnv(path = path, variable = VARS[i], lag = k, year = YR)
          } else {
            getAustraliaEnv(path = path, variable = VARS[i], lag = k, year = DIST[m, 1], month = DIST[m, 2])
          }
        }))
      })
      RAST1 <- do.call(c, VALS[c(1:5, 7, 11)])
      RAST2 <- do.call(c, VALS[c(6, 9:10)])
      RAST3 <- VALS[[8]]
      RAST1 <- crop(RAST1, RAST2)
      RAST3 <- resample(RAST3, RAST1[[1]])
      VALS <- c(RAST1, RAST2, RAST3)
      retn <- data.frame(terra::extract(VALS, pnt))
      colnames(retn) <- c("ID", "temp_max", "tempMax_2mon", "tempMax_3mon",
                          "tempMax_6mon", "tempMax_9mon", "tempMax_12mon",
                          "tempMax_15mon", "tempMax_18mon", "tempMax_21mon",
                          "tempMax_24mon", "temp_min", "tempMin_2mon", 
                          "tempMin_3mon", "tempMin_6mon", "tempMin_9mon",
                          "tempMin_12mon", "tempMin_15mon", "tempMin_18mon",
                          "tempMin_21mon", "tempMin_24mon", "temp_diff", "tempDiff_2mon",
                          "tempDiff_3mon", "tempDiff_6mon", "tempDiff_9mon", 
                          "tempDiff_12mon", "tempDiff_15mon", "tempDiff_18mon",
                          "tempDiff_21mon", "tempDiff_24mon", "prec", "prec_2mon",
                          "prec_3mon", "prec_6mon", "prec_9mon", "prec_12mon",
                          "prec_15mon", "prec_18mon", "prec_21mon", "prec_24mon",
                          "prec_anom", "prec_2anom", "prec_3anom", "prec_6anom",
                          "prec_9anom", "prec_12anom", "prec_15anom", "prec_18anom",
                          "prec_21anom", "prec_24anom", "vaporPressure", "perc_pasture",
                          "perc_urban", "perc_forest", "perc_crops", "ndvi", "ndvi_1mon",
                          "ndvi_2mon", "ndvi_3mon", "ndvi_4mon", "soilMoisture",
                          "evapPotential", "solarExposure", "solarExposure_2mon",
                          "solarExposure_3mon", "solarExposure_6mon",
                          "solarExposure_9mon", "solarExposure_12mon", 
                          "solarExposure_15mon", "solarExposure_18mon",
                          "solarExposure_21mon", "solarExposure_24mon")
      retn$ID <- nums
      if(m == 1) {
        dretn <- retn
      } else dretn <- rbind(dretn, retn)
    }
    dretn <- dretn[order(dretn$ID), ]
    ONI <- read.table(paste0(path, "oni.data.txt"), sep = "", row.names = 1)
    colnames(ONI) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
                       "Sep", "Oct", "Nov", "Dec")
    ONI <- t(as.matrix(ONI))
    ONI[ONI == -99.90] <- NA
    NUMS <- sapply(1:nrow(dretn), function(i) (which(year[i] == 1950:2021) - 1) * 12 + month[i]) 
    ONIm <- sapply(NUMS, function(i) ONI[i])
    ONI_2mon <- sapply(NUMS, function(i) ONI[i - 2])
    ONI_3mon <- sapply(NUMS, function(i) ONI[i - 3])
    ONI_6mon <- sapply(NUMS, function(i) ONI[i - 6])
    ONI_9mon <- sapply(NUMS, function(i) ONI[i - 9])
    ONI_12mon <- sapply(NUMS, function(i) ONI[i - 12])
    ONI_15mon <- sapply(NUMS, function(i) ONI[i - 15])
    ONI_18mon <- sapply(NUMS, function(i) ONI[i - 18])
    ONI_21mon <- sapply(NUMS, function(i) ONI[i - 21])
    ONI_24mon <- sapply(NUMS, function(i) ONI[i - 24])
    dretn <-mutate(dretn, ONI = ONIm,
                   ONI_2mon = ONI_2mon,
                   ONI_3mon = ONI_3mon,
                   ONI_6mon = ONI_6mon,
                   ONI_9mon = ONI_9mon,
                   ONI_12mon = ONI_12mon,
                   ONI_15mon = ONI_15mon,
                   ONI_18mon = ONI_18mon,
                   ONI_21mon = ONI_21mon,
                   ONI_24mon = ONI_24mon)
    SOI <- read.csv(paste0(path, "soi.csv"), row.names = 1)
    SOI <- t(as.matrix(SOI))
    NUMS <- sapply(1:nrow(dretn), function(i) (which(year[i] == 1951:2022) - 1) * 12 + month[i])
    SOIm <- sapply(NUMS, function(i) SOI[i])
    SOI_2mon <- sapply(NUMS, function(i) SOI[i - 2])
    SOI_3mon <- sapply(NUMS, function(i) SOI[i - 3])
    SOI_6mon <- sapply(NUMS, function(i) SOI[i - 6])
    SOI_9mon <- sapply(NUMS, function(i) SOI[i - 9])
    SOI_12mon <- sapply(NUMS, function(i) SOI[i - 12])
    SOI_15mon <- sapply(NUMS, function(i) SOI[i - 15])
    SOI_18mon <- sapply(NUMS, function(i) SOI[i - 18])
    SOI_21mon <- sapply(NUMS, function(i) SOI[i - 21])
    SOI_24mon <- sapply(NUMS, function(i) SOI[i - 24])
    dretn <-mutate(dretn, SOI = SOIm,
                   SOI_2mon = SOI_2mon,
                   SOI_3mon = SOI_3mon,
                   SOI_6mon = SOI_6mon,
                   SOI_9mon = SOI_9mon,
                   SOI_12mon = SOI_12mon,
                   SOI_15mon = SOI_15mon,
                   SOI_18mon = SOI_18mon,
                   SOI_21mon = SOI_21mon,
                   SOI_24mon = SOI_24mon)
    SAM <- read.table(paste0(path, "newsam.1957.2007.txt"), fill = T)
    SAM <- t(as.matrix(SAM))
    NUMS <- sapply(1:nrow(dretn), function(i) (which(year[i] == 1957:2022) - 1) * 12 + month[i])
    SAMm <- sapply(NUMS, function(i) SAM[i])
    SAM_2mon <- sapply(NUMS, function(i) SAM[i - 2])
    SAM_3mon <- sapply(NUMS, function(i) SAM[i - 3])
    SAM_6mon <- sapply(NUMS, function(i) SAM[i - 6])
    SAM_9mon <- sapply(NUMS, function(i) SAM[i - 9])
    SAM_12mon <- sapply(NUMS, function(i) SAM[i - 12])
    SAM_15mon <- sapply(NUMS, function(i) SAM[i - 15])
    SAM_18mon <- sapply(NUMS, function(i) SAM[i - 18])
    SAM_21mon <- sapply(NUMS, function(i) SAM[i - 21])
    SAM_24mon <- sapply(NUMS, function(i) SAM[i - 24])
    dretn <-mutate(dretn, SAM = SAMm,
                   SAM_2mon = SAM_2mon,
                   SAM_3mon = SAM_3mon,
                   SAM_6mon = SAM_6mon,
                   SAM_9mon = SAM_9mon,
                   SAM_12mon = SAM_12mon,
                   SAM_15mon = SAM_15mon,
                   SAM_18mon = SAM_18mon,
                   SAM_21mon = SAM_21mon,
                   SAM_24mon = SAM_24mon)
    return(dretn)
  }
}
