#Model of number of BFF roosts per quarter from 2008 - 2021

#Generate random draws from normal distribution centered on fitted value with sd=standard error
#Truncate lower bound at 1
#Write as function to return matrix with columns equal to the number of replicates needed
roost.counts <- function(n, months="quarter"){ 
  #n = number of random replicates to generate for the time period above (2008-2021/1)
  # months = "quarter" gives the same value for all 3 months in the quarter
  bff.roosts <- read.csv("https://raw.githubusercontent.com/hanlab-ecol/BatOneHealth/main/data/bff_one_week_per_quarter_illya.csv")
  bff.roosts <- cbind(bff.roosts,as.numeric(substring(bff.roosts$year_quarter,1,4)),
                      as.numeric(substring(bff.roosts$year_quarter,6,6)),
                      as.numeric(as.factor(bff.roosts$year_quarter)))
  names(bff.roosts)[6:8] <- c("Year","Quarter","Time")
  roosts.num <- bff.roosts[7:(dim(bff.roosts)[1]-1),] #Take off incomplete years at beginning and zero value at end
  #Starts Jan 2008
  #Following: https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/
  ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
  ## No AR based on model comparison
  roost.gam <- mgcv::gamm(bff.presence.count ~ s(Quarter, bs = "cc", k = 4) + s(Time, k = 53),
                          data = roosts.num,
                          control = ctrl)
  #This version produces a different roost number FOR EACH MONTH
  if (months!="quarter"){
    Time <- rep(7:59,each=3) #These are the roost values monthly with the same values for each month in that quarter
    Quarter <- rep(rep(c(1:4),each=3),length(unique(roosts.num$Year))) #Quarter value for each year
    Quarter <- Quarter[1:length(Time)] #Trim extra values for incomplete years at the end of the data
    pred.points <- data.frame(Quarter,Time)
    pred.roosts <- predict(roost.gam$gam,  newdata = pred.points, se.fit = TRUE)

    n.roosts <- replicate(n, rtruncnorm(n=length(Time), a=1, mean=pred.roosts$fit, sd=pred.roosts$se.fit))
    #return(n.roosts)
  } else { #This repeats the same value for each month in a quarter
    Time <- 7:59
    Quarter <- rep(c(1:4),ceiling(length(Time)/4))
    Quarter <- Quarter[1:length(Time)] #Trim extra values for incomplete years at the end of the data
    pred.points <- data.frame(Quarter,Time)
    pred.roosts <- predict(roost.gam$gam,  newdata = pred.points, se.fit = TRUE)
    n.roosts <- replicate(n, truncnorm::rtruncnorm(n=length(Time), a=1, mean=pred.roosts$fit, sd=pred.roosts$se.fit))
    n.roosts <- apply(X=n.roosts, MARGIN=2, FUN=rep, each=3)
    #return(n.roosts)
  }
  anchor <- cbind(year = rep(dplyr::distinct(roosts.num, Year, Quarter)$Year, each = 3),
                   month = rep(1:12, length(unique(roosts.num$Year)))[1:159])
  return(cbind(anchor, n.roosts))
  #This version produces the same roost number for all three months in the quarter
}

# plot1 <- ggplot(roosts.num, aes(x=Time, y=bff.presence.count)) +
#   geom_point() +
#   geom_line(aes(y=predict)) +
#   xlab("Time by quarter 2008-2021") +
#   #scale_x_discrete(breaks=roosts.num$Time,label=roosts.num$year_quarter) +
#   theme_minimal() +
#   theme(axis.text.x=element_text(angle=80, hjust=1)) 
