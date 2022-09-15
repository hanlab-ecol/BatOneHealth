#Load required packages
library(raster)
library(gbm)
library(sf)
library(terra)
library(dplyr)
library(doRNG)
library(ggplot2)
library(mgcv)
library(truncnorm)
setwd("~/Documents/Bat spillover") #Adjust for computer
#Functions for accessing and attaching environmental data
source("https://raw.githubusercontent.com/hanlab-ecol/BatOneHealth/main/getAustraliaEnv.R") #Used in next function to access specific data frames
source("https://raw.githubusercontent.com/hanlab-ecol/BatOneHealth/main/attachAustraliaEnv.R") #Used to get the needed data to run various models
source("https://raw.githubusercontent.com/hanlab-ecol/BatOneHealth/main/estimatePrevalence.R") #Used to get output, requires above functions and accessory functions below
#All raster directories need to be in the working directory
source("https://raw.githubusercontent.com/hanlab-ecol/BatOneHealth/main/roost.counts")

#NECESSARY FUNCTIONS
#Select roost locations, function called later 
r.location <- function(suitability, n){
  roosts <- spatSample(suitability,n,method="weights",
             replace=FALSE,na.rm=TRUE,as.points=TRUE)#Returns location(s)
}

#function to estimate rehab probability
rehab.prob <- function(roost.data, rehab.model){
  DATA <- as.data.frame(roost.data)
  rehab <- predict(rehab.model, newdata=DATA, type="response", verbose=FALSE)
  roost.data <- cbind(roost.data,rehab)
  return(roost.data)
}

#function to estimate new roost probability
new.roost.prob <- function(roost.data, new.roost.model){
  DATA <- as.data.frame(roost.data)
  new.roost <- predict(new.roost.model,newdata=DATA, type="response",verbose=FALSE)
  roost.data <- cbind(roost.data,new.roost)
  return(roost.data)
}

food.shortage.prob <- function(food.shortage.pred, date, lag_months){
  date <- date
  num <- which(food.shortage.pred$Date==date)
  return(mean(food.shortage.pred[(num - lag_months):num, "Probability"]))
}

#output only returned in memory if return.stack==TRUE
cl <- parallel::makeCluster(7, "SOCK") #used to run parallel processing, ignore if not wanted
doSNOW::registerDoSNOW(cl) #same as above line
output <-estimatePrevalence(years=2008:2019,
                            months=1:12,
                            reps=100,
                            buffer = 50000,
                            ref.prevalence=0.66,
                            lag_months = 12,
                            models = c("food shortage", "rehab", "new roost"),
                            averaging = NA,
                            cl = cl,
                            return.stack=FALSE)
#The above will write raster data to your working directory
