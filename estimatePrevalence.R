
#'This function operates to generate predicted Hendra prevalence in the form of monthly
#' raster layers.
#' @param path This is the path to a repository with all of the environmental data and fit models used in this project. By default, the function will grab the current working directory.
#' @param years The years for which to estimate prevalence. Default is 2008-2019. This is full data for food shortage in current table and those months that we have full roost count information for.
#' @param months The months for which to estimate prevalence. Default is 1-12.
#' @param reps The number of replicate, stochastic draws of roost locations for a given month. Default is 100.
#' @param buffer The maximum size allowed for foraging areas associated with individual roosts. Default is radius = 50km. This = 50000m.
#' @param ref.prevalence The maximum prevalence given stressed bats. The default is 0.66 from observed data.
#' @param lag_months The length of the lag to include in calculations of food shortage probability, creating an average of food stress over a period of time. The default is 12 months and a value of 0 months will make sure that no average is used.
#' @param models Dictates the component models to use for estimating prevalence in this multi scale model. By default all three models (food shortage, new roost, and rehab) are included.
#' @param averaging An argument that allows one to essentially downweight rarer areas not chosen as often. The default value is NA, but switching to 0 will reclassify NAs in the raster to 0, giving rarely chosen areas lower prevalence scores when averaging across repetitions
#' @param return.stack Whether or not to retain all the annual raster stacks in memory
#' @param cl A cluster object that can be provided to parallelize this function across the years. The default is NULL, which will run everything on a single core.
#' @return This function will return an object of terra::SpatRaster-class ranging from 
#' 1 to 12 layers representing the variable chosen at the time points (year and 
#' possibly month) chosen. Default is to save to file on 1 year timesteps.
#' @example estimatePrevalence(path = "~/Desktop/AustralianClimateData", n.roosts=104, years = c(1996:2019), months = c(1:12), buffer = 50000, prevalence=0.24)

estimatePrevalence <- function(path=getwd(), years=c(2008:2019), months=c(1:12), reps= 100, buffer=50000, ref.prevalence = 0.66, lag_months = 12, models = c("food shortage", "new roost", "rehab"), averaging = NA, return.stack=FALSE, cl = NULL){
  if(sum(!models %in% c("food shortage", "new roost", "rehab")) > 0) stop("One or more models specified are not compatible with this function. Please check spelling and the documentation for accepted models.")
  #Read in the model for bat rehab probability
  load("RehabModel_01062022.Rdata") #REHB
  #Read in new overwintering model
  load("NewRoostPredictedSurfacesModel.Rdata") #new.roost.model
  #Read in food shortage data: this currently goes to 2020-02
  food.shortage.pred <- read.csv("shortage_reduced_environmental.txt")
  names(food.shortage.pred) <- c("Date","Probability")
  set.seed=5302022
  n.roosts <- round(roost.counts(n = reps, months = "quarter"), 0)
  #Loop through years and months
  if(!is.null(cl)) {
  out <- foreach(i = 1:length(years), .packages = c("terra", "sf", "dplyr", "gbm"), .export = c("r.location", "rehab.prob", "new.roost.prob", "getAustraliaEnv", "attachAustraliaEnv", "food.shortage.prob")) %dorng% {
    #Get the SDMs for the given year
    roosts.rast <- rast(paste("AustralianClimateData/BatRoostPredictions/BatRoostPredictions_",years[i],".tif",sep=""))   #SpatRaster for randomization
    #roosts.rast <- project(roosts.rast,y="epsg:8059") #Tried converting raster instead of points below but didn't work
    #Set up template raster for inputting prevalence later
    template <- rast(roosts.rast[[1]])
    #template[!is.na(template)] <- 0 #I think unnecessary, part of older version
    out <- template
    year.stack <- template
    for (j in months){
      #First get the SDM for roosts and set NA values to 0
      roost.presence <- roosts.rast[[j]]
      roost.presence[is.nan(roost.presence)] <- 0
      n.roost <- n.roosts[n.roosts[, 1] == years[i] & n.roosts[, 2] == j, -1:-2]
      #n.roosts <- ifelse(length(n.roosts)>1,n.roosts[((i-1)*12)+j],n.roosts) #This takes the number of roosts for the month if using, otherwise the fix number
      month <- ifelse(j<10,paste(0,j,sep=""),j)
      time <- paste(years[i],"-",month,sep="")
      
      #Calculate a single replicate for roosts, host condition at roosts, and prevalence
      #Accumulate replicates in a raster stack
      stack <- template
      path=paste(getwd(),"/AustralianClimateData/",sep="") #This will be the path that contains the environmental data folders
      # for (k in 1:reps){
      roosts <- lapply(n.roost, function(k) r.location(roost.presence, n = k))
      #roosts <- replicate(reps, r.location(roost.presence,n=n.roosts)) #Uses helper function to randomize roost location weighted by the SDM
      #Get environmental data for the roosts
      roosts1 <- lapply(roosts, function(k) sf::st_as_sf(k)) #This is in sf format, used in this projection and format below
      roosts2 <- lapply(roosts1, function(k) sf::st_transform(k, crs=8059)) #This gets roosts in the same projection as the environmental data (Unsure why formats differ initially, seems subtle). Used only in extracting data.
      RNGE <- c(1:12, 1:12) #creates range of months
      RNUM <- which(RNGE == j)[2] #grabs second appearance of the month in question to easily query across years
      mnths <- RNGE[(RNUM - 11):RNUM] #grabs months in question for lag
      yrs <- c(rep(years[i] - 1, 12 - j), rep(years[i], j)) #grabs years to match lag months
      roosts3 <- lapply(roosts2, function(k) rep(list(k), 12)) #replicates roost data 12 times
      roosts3 <- lapply(roosts3, function(k) do.call(rbind, lapply(1:12, function(p) mutate(k[[p]], month = mnths[p], year = yrs[p]))))
      roost.data <- attachAustraliaEnv(points = do.call(rbind, roosts3), path=path, year = do.call(rbind, roosts3)$year, month = do.call(rbind, roosts3)$month, log.transform = c("tempMax_2mon", "tempMax_3mon", "prec", "prec_2mon", "tempDiff_2mon","prec_2anom","prec_3anom","solarExposure","solarExposure_12mon"))
      roost.data <- cbind(roost.data, index = do.call(c, sapply(1:reps, function(k) rep(k, each = n.roost[k] * 12))), row_num = do.call(c, sapply(1:reps, function(p) rep(1:n.roost[p], 12)))) #add the rep number to each point according to the number of roosts #multiply by lag months to get the appropriate index
      roost.data <- lapply(1:reps, function(k) roost.data[roost.data$index == k, ]) #return to list format
      #probability of bat rehab at roosts, column appended
      if(sum(models %in% "rehab") == 1) {
        roost.data <- lapply(roost.data, function(k) rehab.prob(roost.data = k, rehab.model = REHB))
      }
      #probability of new roost state at roosts, column appended
      if(sum(models %in% "new roost") == 1) {
        roost.data <- lapply(roost.data, function(k) new.roost.prob(roost.data = k, new.roost.model = SURF))
      }
      #food stress condition for the month in question. Data is in YYYY-MM format, formatted above
      if(sum(models %in% "food shortage") == 1) {
        food.short <- food.shortage.prob(food.shortage.pred=food.shortage.pred, date=time, lag_months = lag_months)
      }
      #Combine stress factors into an index. The bat rehab and new roosts will be averaged and then combined additively with the food stress
      if(sum(models %in% c("rehab", "new roost")) == 2) {
        roost.data <- lapply(roost.data, function(k) {
          k %>%
            group_by(row_num) %>%
            summarise(rehab = mean(rehab, na.rm = T),
                      new.roost = mean(new.roost, na.rm = T)) %>%
            ungroup() %>%
            data.frame()
        })
        space.stress <- lapply(roost.data, function(k) sqrt(k$rehab * k$new.roost))
      } else if (sum(models %in% c("rehab", "new roost")) == 1 & sum(models %in% "rehab") == 1) {
        roost.data <- lapply(roost.data, function(k) {
          k %>%
            group_by(row_num) %>%
            summarise(rehab = mean(rehab, na.rm = T)) %>%
            ungroup() %>%
            data.frame()
        })
        space.stress <- lapply(roost.data, function(k) k$rehab)
      } else if (sum(models %in% c("rehab", "new roost")) == 1 & sum(models %in% "new roost") == 1) {
        roost.data <- lapply(roost.data, function(k) {
          k %>%
            group_by(row_num) %>%
            summarise(new.roost = mean(new.roost, na.rm = T)) %>%
            ungroup() %>%
            data.frame()
        })
        space.stress <- lapply(roost.data, function(k) k$new.roost)
      }
      if(sum(models %in% c("rehab", "new roost")) > 0 & sum(models %in% "food shortage") == 1) {
        stress <- lapply(space.stress, function(k) 1 - ((1 - k) * (1 - food.short)))
      } else if (sum(models %in% c("rehab", "new roost")) > 0 & sum(models %in% "food shortage") == 0) {
        stress <- space.stress
      } else if (sum(models %in% c("rehab", "new roost")) == 0) {
        stress <- lapply(n.roost, function(k) rep(food.short, k))
      }
      #Multiply index by reference prevalence, this represents maximum prevalence when stress is 100% likely
      prevalence <- lapply(stress, function(k) k * ref.prevalence)
      #roost.data2 <- cbind(roost.data,food.short,space.stress,stress,prevalence) #Currently this is not used or output
      #Place prevalence on map by doing tesselation and then clipping by foraging area size
      v <- lapply(roosts, function(k) terra::voronoi(k))
      v1 <- lapply(v, function(k) sf::st_as_sf(k))
      #Need to do individual unions, but the tessellation and buffers are not in the same order, get intersects first to make sure the polygons match the points
      inters <- lapply(1:length(roosts1), function(k) unlist(sf::st_intersects(v1[[k]], roosts1[[k]])))
      v1 <- lapply(1:length(v1), function(k) cbind(v1[[k]], index = inters[[k]]))
      buffer.50k<- lapply(roosts1, function(k) sf::st_buffer(k, dist = buffer)) #buffer each roost by 50km]
      buffer.50k2 <- lapply(1:reps, function(k) buffer.50k[[k]][inters[[k]], ])
      buffer.50k2 <- lapply(1:length(buffer.50k2), function(k) cbind(buffer.50k2[[k]], index = inters[[k]])) #This index matches v1 to facilitate matching intersected polygons below
      for(k in v1) st_agr(k) = "constant"
      for(k in buffer.50k2) st_agr(k) = "constant"
      map <- lapply(1:length(v1), function(k) st_intersection(v1[[k]], buffer.50k2[[k]]))
      map <- lapply(map, function(k) k[which(k[, 2] == k[ ,4]), ])
      map <- lapply(map, function(k) k[!st_is_empty(k), ])
      if(sum(models %in% c("rehab", "new roost")) == 2) {
        map <- lapply(1:reps, function(k) cbind(map[[k]], prevalence = prevalence[[k]][inters[[k]]], rehab = roost.data[[k]]$rehab[inters[[k]]], newroosts = roost.data[[k]]$new.roost[inters[[k]]])) #This is the prevalence map for iteration k of the roost locations, order is to match roost locations
      } else if (sum(models %in% c("rehab", "new roost")) == 1 & sum(models %in% "rehab") == 1) {
        map <- lapply(1:reps, function(k) cbind(map[[k]], prevalence = prevalence[[k]][inters[[k]]], rehab = roost.data[[k]]$rehab[inters[[k]]]))
      } else if (sum(models %in% c("rehab", "new roost")) == 1 & sum(models %in% "new roost") == 1) {
        map <- lapply(1:reps, function(k) cbind(map[[k]], prevalence = prevalence[[k]][inters[[k]]], newroosts = roost.data[[k]]$new.roost[inters[[k]]]))
      } else {
        map <- lapply(1:reps, function(k) cbind(map[[k]], prevalence = prevalence[[k]][inters[[k]]]))
      }
      #Transfer values to raster
      rast.map <- lapply(1:reps, function(k) rasterize(vect(map[[k]]), template, "prevalence"))
      if(!is.na(averaging)) rast.map <- lapply(rast.map, function(k) classify(k, cbind(NA, 0)))
      rast.map <- lapply(rast.map, function(k) mask(k,roosts.rast[[1]])) #Takes out areas that are NA in the underlying environmental data: i.e. restricts to coastline
      if(sum(models %in% "rehab") == 1) {
        rehab.map <- lapply(1:reps, function(k) rasterize(vect(map[[k]]), template, "rehab"))
        if(!is.na(averaging)) rehab.map <- lapply(rehab.map, function(k) classify(k, cbind(NA, 0)))
        rehab.map <- lapply(rehab.map, function(k) mask(k, roosts.rast[[1]]))
      }
      if(sum(models %in% "new roost") == 1) {
        newroosts.map <- lapply(1:reps, function(k) rasterize(vect(map[[k]]), template, "newroosts"))
        if(!is.na(averaging)) newroosts.map <- lapply(newroosts.map, function(k) classify(k, cbind(NA, 0)))
        newroosts.map <- lapply(newroosts.map, function(k) mask(k, roosts.rast[[1]]))
      }
      stack <- do.call(c, rast.map) #stack the k rasters
      if(sum(models %in% "rehab") == 1) rehab.stack <- do.call(c, rehab.map)
      if(sum(models %in% "new roost") == 1) newroosts.stack <- do.call(c, newroosts.map)
      message("Completed",years[i],"month",j)
      year.stack <- c(year.stack,terra::app(stack,mean,na.rm=T)) #This averages the k rasters for the month and adds as a layer to the annual stack
      if(j == 1) {
        if(sum(models %in% "rehab") == 1) year.rehab <- terra::app(rehab.stack, mean, na.rm = T)
        if(sum(models %in% "new roost") == 1) year.newroost <- terra::app(newroosts.stack, mean, na.rm = T)
      }  else {
        if(sum(models %in% "rehab") == 1) year.rehab <- c(year.rehab, terra::app(rehab.stack, mean, na.rm = T))
        if(sum(models %in% "new roost") == 1) year.newroost <- c(year.newroost, terra::app(newroosts.stack, mean, na.rm = T))
      }
      # }
    }
    #Write the annual stack to a file
    writeRaster(year.stack,filename=paste("ExpectedPrevalence_",years[i],".tif",sep=""),overwrite=TRUE)
    if(sum(models %in% "rehab") == 1) writeRaster(year.rehab, filename = paste0("RehabProbability_", years[i], ".tif"), overwrite = T)
    if(sum(models %in% "new roost") == 1) writeRaster(year.newroost, filename = paste0("NewRoostProbability_", years[i], ".tif"), overwrite = T)
    if (return.stack==TRUE) c(out,year.stack)
  }
  }
  else {
    for(i in 1:length(years)) {
      #Get the SDMs for the given year
      roosts.rast <- rast(paste("AustralianClimateData/BatRoostPredictions/BatRoostPredictions_",years[i],".tif",sep=""))   #SpatRaster for randomization
      #roosts.rast <- project(roosts.rast,y="epsg:8059") #Tried converting raster instead of points below but didn't work
      #Set up template raster for inputting prevalence later
      template <- rast(roosts.rast[[1]])
      #template[!is.na(template)] <- 0 #I think unnecessary, part of older version
      out <- template
      year.stack <- template
      for (j in months){
        #First get the SDM for roosts and set NA values to 0
        roost.presence <- roosts.rast[[j]]
        roost.presence[is.nan(roost.presence)] <- 0
        n.roost <- n.roosts[n.roosts[, 1] == years[i] & n.roosts[, 2] == j, -1:-2]
        #n.roosts <- ifelse(length(n.roosts)>1,n.roosts[((i-1)*12)+j],n.roosts) #This takes the number of roosts for the month if using, otherwise the fix number
        month <- ifelse(j<10,paste(0,j,sep=""),j)
        time <- paste(years[i],"-",month,sep="")
        
        #Calculate a single replicate for roosts, host condition at roosts, and prevalence
        #Accumulate replicates in a raster stack
        stack <- template
        path=paste(getwd(),"/AustralianClimateData/",sep="") #This will be the path that contains the environmental data folders
        # for (k in 1:reps){
        roosts <- lapply(n.roost, function(k) r.location(roost.presence, n = k))
        #roosts <- replicate(reps, r.location(roost.presence,n=n.roosts)) #Uses helper function to randomize roost location weighted by the SDM
        #Get environmental data for the roosts
        roosts1 <- lapply(roosts, function(k) sf::st_as_sf(k)) #This is in sf format, used in this projection and format below
        roosts2 <- lapply(roosts1, function(k) sf::st_transform(k, crs=8059)) #This gets roosts in the same projection as the environmental data (Unsure why formats differ initially, seems subtle). Used only in extracting data.
        RNGE <- c(1:12, 1:12) #creates range of months
        RNUM <- which(RNGE == j)[2] #grabs second appearance of the month in question to easily query across years
        mnths <- RNGE[(RNUM - 11):RNUM] #grabs months in question for lag
        yrs <- c(rep(years[i] - 1, 12 - j), rep(years[i], j)) #grabs years to match lag months
        roosts3 <- lapply(roosts2, function(k) rep(list(k), 12)) #replicates roost data 12 times
        roosts3 <- lapply(roosts3, function(k) do.call(rbind, lapply(1:12, function(p) mutate(k[[p]], month = mnths[p], year = yrs[p]))))
        roost.data <- attachAustraliaEnv(points = do.call(rbind, roosts3), path=path, year = do.call(rbind, roosts3)$year, month = do.call(rbind, roosts3)$month, log.transform = c("tempMax_2mon", "tempMax_3mon", "prec", "prec_2mon", "tempDiff_2mon","prec_2anom","prec_3anom","solarExposure","solarExposure_12mon"))
        roost.data <- cbind(roost.data, index = do.call(c, sapply(1:reps, function(k) rep(k, each = n.roost[k] * 12))), row_num = do.call(c, sapply(1:reps, function(p) rep(1:n.roost[p], 12)))) #add the rep number to each point according to the number of roosts
        roost.data <- lapply(1:reps, function(k) roost.data[roost.data$index == k, ]) #return to list format
        #probability of bat rehab at roosts, column appended
        if(sum(models %in% "rehab") == 1) {
          roost.data <- lapply(roost.data, function(k) rehab.prob(roost.data = k, rehab.model = REHB))
        }
        #probability of new roost state at roosts, column appended
        if(sum(models %in% "new roost") == 1) {
          roost.data <- lapply(roost.data, function(k) new.roost.prob(roost.data = k, new.roost.model = SURF))
        }
        #food stress condition for the month in question. Data is in YYYY-MM format, formatted above
        if(sum(models %in% "food shortage") == 1) {
          food.short <- food.shortage.prob(food.shortage.pred=food.shortage.pred, date=time, lag_months = lag_months)
        }
        #Combine stress factors into an index. The bat rehab and new roosts will be averaged and then combined additively with the food stress
        if(sum(models %in% c("rehab", "new roost")) == 2) {
          roost.data <- lapply(roost.data, function(k) {
            k %>%
              group_by(row_num) %>%
              summarise(rehab = mean(rehab, na.rm = T),
                        new.roost = mean(new.roost, na.rm = T)) %>%
              ungroup() %>%
              data.frame()
          })
          space.stress <- lapply(roost.data, function(k) sqrt(k$rehab * k$new.roost))
        } else if (sum(models %in% c("rehab", "new roost")) == 1 & sum(models %in% "rehab") == 1) {
          roost.data <- lapply(roost.data, function(k) {
            k %>%
              group_by(row_num) %>%
              summarise(rehab = mean(rehab, na.rm = T)) %>%
              ungroup() %>%
              data.frame()
          })
          space.stress <- lapply(roost.data, function(k) k$rehab)
        } else if (sum(models %in% c("rehab", "new roost")) == 1 & sum(models %in% "new roost") == 1) {
          roost.data <- lapply(roost.data, function(k) {
            k %>%
              group_by(row_num) %>%
              summarise(new.roost = mean(new.roost, na.rm = T)) %>%
              ungroup() %>%
              data.frame()
          })
          space.stress <- lapply(roost.data, function(k) k$new.roost)
        }
        if(sum(models %in% c("rehab", "new roost")) > 0 & sum(models %in% "food shortage") == 1) {
          stress <- lapply(space.stress, function(k) 1 - ((1 - k) * (1 - food.short)))
        } else if (sum(models %in% c("rehab", "new roost")) > 0 & sum(models %in% "food shortage") == 0) {
          stress <- space.stress
        } else if (sum(models %in% c("rehab", "new roost")) == 0) {
          stress <- lapply(n.roost, function(k) rep(food.short, k))
        }
        #Multiply index by reference prevalence, this represents maximum prevalence when stress is 100% likely
        prevalence <- lapply(stress, function(k) k * ref.prevalence)
        #roost.data2 <- cbind(roost.data,food.short,space.stress,stress,prevalence) #Currently this is not used or output
        #Place prevalence on map by doing tesselation and then clipping by foraging area size
        v <- lapply(roosts, function(k) terra::voronoi(k))
        v1 <- lapply(v, function(k) sf::st_as_sf(k))
        #Need to do individual unions, but the tessellation and buffers are not in the same order, get intersects first to make sure the polygons match the points
        inters <- lapply(1:length(roosts1), function(k) unlist(sf::st_intersects(v1[[k]], roosts1[[k]])))
        v1 <- lapply(1:length(v1), function(k) cbind(v1[[k]], index = inters[[k]]))
        buffer.50k<- lapply(roosts1, function(k) sf::st_buffer(k, dist = buffer)) #buffer each roost by 50km]
        buffer.50k2 <- lapply(1:reps, function(k) buffer.50k[[k]][inters[[k]], ])
        buffer.50k2 <- lapply(1:length(buffer.50k2), function(k) cbind(buffer.50k2[[k]], index = inters[[k]])) #This index matches v1 to facilitate matching intersected polygons below
        for(k in v1) st_agr(k) = "constant"
        for(k in buffer.50k2) st_agr(k) = "constant"
        map <- lapply(1:length(v1), function(k) st_intersection(v1[[k]], buffer.50k2[[k]]))
        map <- lapply(map, function(k) k[which(k[, 2] == k[ ,4]), ])
        map <- lapply(map, function(k) k[!st_is_empty(k), ])
        if(sum(models %in% c("rehab", "new roost")) == 2) {
          map <- lapply(1:reps, function(k) cbind(map[[k]], prevalence = prevalence[[k]][inters[[k]]], rehab = roost.data[[k]]$rehab[inters[[k]]], newroosts = roost.data[[k]]$new.roost[inters[[k]]])) #This is the prevalence map for iteration k of the roost locations, order is to match roost locations
        } else if (sum(models %in% c("rehab", "new roost")) == 1 & sum(models %in% "rehab") == 1) {
          map <- lapply(1:reps, function(k) cbind(map[[k]], prevalence = prevalence[[k]][inters[[k]]], rehab = roost.data[[k]]$rehab[inters[[k]]]))
        } else if (sum(models %in% c("rehab", "new roost")) == 1 & sum(models %in% "new roost") == 1) {
          map <- lapply(1:reps, function(k) cbind(map[[k]], prevalence = prevalence[[k]][inters[[k]]], newroosts = roost.data[[k]]$new.roost[inters[[k]]]))
        } else {
          map <- lapply(1:reps, function(k) cbind(map[[k]], prevalence = prevalence[[k]][inters[[k]]]))
        }
        #Transfer values to raster
        rast.map <- lapply(1:reps, function(k) rasterize(vect(map[[k]]), template, "prevalence"))
        if(!is.na(averaging)) rast.map <- lapply(rast.map, function(k) classify(k, cbind(NA, 0)))
        rast.map <- lapply(rast.map, function(k) mask(k,roosts.rast[[1]])) #Takes out areas that are NA in the underlying environmental data: i.e. restricts to coastline
        if(sum(models %in% "rehab") == 1) {
          rehab.map <- lapply(1:reps, function(k) rasterize(vect(map[[k]]), template, "rehab"))
          if(!is.na(averaging)) rehab.map <- lapply(rehab.map, function(k) classify(k, cbind(NA, 0)))
          rehab.map <- lapply(rehab.map, function(k) mask(k, roosts.rast[[1]]))
        }
        if(sum(models %in% "new roost") == 1) {
          newroosts.map <- lapply(1:reps, function(k) rasterize(vect(map[[k]]), template, "newroosts"))
          if(!is.na(averaging)) newroosts.map <- lapply(newroosts.map, function(k) classify(k, cbind(NA, 0)))
          newroosts.map <- lapply(newroosts.map, function(k) mask(k, roosts.rast[[1]]))
        }
        stack <- do.call(c, rast.map) #stack the k rasters
        if(sum(models %in% "rehab") == 1) rehab.stack <- do.call(c, rehab.map)
        if(sum(models %in% "new roost") == 1) newroosts.stack <- do.call(c, newroosts.map)
        message("Completed",years[i],"month",j)
        year.stack <- c(year.stack,terra::app(stack,mean,na.rm=T)) #This averages the k rasters for the month and adds as a layer to the annual stack
        if(j == 1) {
          if(sum(models %in% "rehab") == 1) year.rehab <- terra::app(rehab.stack, mean, na.rm = T)
          if(sum(models %in% "new roost") == 1) year.newroost <- terra::app(newroosts.stack, mean, na.rm = T)
        }  else {
          if(sum(models %in% "rehab") == 1) year.rehab <- c(year.rehab, terra::app(rehab.stack, mean, na.rm = T))
          if(sum(models %in% "new roost") == 1) year.newroost <- c(year.newroost, terra::app(newroosts.stack, mean, na.rm = T))
        }
        # }
      }
      #Write the annual stack to a file
      writeRaster(year.stack,filename=paste("ExpectedPrevalence_",years[i],".tif",sep=""),overwrite=TRUE)
      if(sum(models %in% "rehab") == 1) writeRaster(year.rehab, filename = paste0("RehabProbability_", years[i], ".tif"), overwrite = T)
      if(sum(models %in% "new roost") == 1) writeRaster(year.newroost, filename = paste0("NewRoostProbability_", years[i], ".tif"), overwrite = T)
      if (return.stack==TRUE) out <- c(out,year.stack)
    }
  }
  return(out)
}
