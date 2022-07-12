##############
# Comparing prevalence predictions to shedding data
##############
# libraries
library(terra)
library(dplyr)
library(sf)
library(wesanderson)
library(ggplot2)
library(viridis)
library(ggspatial)

################
# importing Hume Fields data
hev_raw = read.csv('data/hendra-virus-transmission-data-east-australia.csv', header = TRUE)
head(hev_raw)
hev_raw$date = as.Date(hev_raw$date, format = '%d/%m/%Y')
hev_prev = as.data.frame(hev_raw %>%
                           group_by(date,loc) %>%
                           summarize(prev = mean(res, na.rm = TRUE),
                                     obs = length(res),
                                     fftotal = mean(fftotal, na.rm = TRUE),
                                     estb = mean(estb, na.rm = TRUE),
                                     pct_b = mean(pct_b, na.rm = TRUE)))

hev_loc = read.csv('data/locations_sites.csv', header = TRUE)
hev_prev_loc = merge(hev_prev, hev_loc, by.x = 'loc', by.y ='site')
range(hev_prev_loc$UTM_zone)
range(hev_prev_loc$northing)

# Two different UTM zones > convert and match
hev_prev_56 = hev_prev_loc[hev_prev_loc$UTM_zone==56,]
hev_prev_55 = hev_prev_loc[hev_prev_loc$UTM_zone==55,]

unique(hev_prev_loc$loc)[order(unique(hev_prev_loc$loc))]
hendra_sf_56 = st_as_sf(hev_prev_56, coords = c( 'easting','northing'),
                        crs = "+proj=utm +zone=56 +south +datum=WGS84")
crs_pred ='+proj=lcc +lat_0=-32 +lon_0=135 +lat_1=-28 +lat_2=-36 +x_0=1000000 +y_0=2000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
hendra_56_crs = st_transform(hendra_sf_56, crs_pred)

hendra_sf_55 = st_as_sf(hev_prev_55, coords = c( 'easting','northing'),
                        crs = "+proj=utm +zone=55 +south +datum=WGS84")
hendra_55_crs = st_transform(hendra_sf_55, crs_pred)
hendra_crs = rbind(hendra_56_crs,hendra_55_crs)


##########
# australia shapefile (can use any)
ausdir = "~/Dropbox/Projects/Hendra/workspace/sp_dat/AUS_adm/nsaasr9nnd_02211a04es_geo___"
aus = st_read(ausdir,'aust_cd66states', crs = "+proj=longlat +datum=WGS84") #country boundary
aus_crs = st_transform(aus, crs_pred)
aus_simp <- st_simplify(aus_crs, preserveTopology = FALSE, 
                        dTolerance = 10000)
# plotting 26 sites
ggplot() + 
  geom_sf(data = aus_crs, lwd = 0.5, col = 'black',fill = 'grey90')+
  geom_sf(data = hendra_crs, aes(col = loc), size = 2) +
  geom_sf_text(data = hendra_crs, aes(label = loc), colour = "black", cex = 1.2)+
  #scale_color_viridis(discrete=TRUE, option = 'B')+
  #scale_color_manual(values = roost_col)+
  coord_sf(xlim = c(2066633, 2903360), ylim = c(1506948, 3977938))+
  theme_bw()+
  annotation_scale(location = "bl", width_hint = 0.2) +
  theme(legend.position="none")


###############
# importing prevalence rasters and transforming

# need to set working directory to prevalence rasters
prev_temp <- list.files("prev_with0", full.names = TRUE)
prev_field <- rast(prev_temp[5:8]) # years overlapping with shedding data

monthly_estimates = as.data.frame(matrix(NA, nrow = nrow(hendra_crs)*nlyr(prev_field), ncol =1))
monthly_estimates$site = rep(hendra_crs$loc, each=nlyr(prev_field))
colnames(monthly_estimates) = c('prevalence_est', 'site')
for(j in 1:nrow(hendra_crs)){
  #j = 1 
  buff =  vect(st_buffer(hendra_crs[j,], dist = 1000))
  for(i in 1:nlyr(prev_field)){
    #i = 10
    values_buff = terra::extract(prev_field[[i]], buff)$mean
    monthly_estimates$prevalence_est[(j-1)*nlyr(prev_field)+i] = mean(values_buff[!is.na(values_buff)]) 
  }
}
monthly_estimates$month = rep(1:12, times=4*nrow(hendra_crs))
monthly_estimates$year = rep(rep(2011:2014, each=12), times =nrow(hendra_crs))
monthly_estimates$day = 15
monthly_estimates$date = as.Date(paste(monthly_estimates$year, monthly_estimates$month, monthly_estimates$day, sep="-"), "%Y-%m-%d")
#monthly_estimates$prevalence_adj= monthly_estimates$prevalence_est*1.67
monthly_estimates$site = as.factor(monthly_estimates$site)
# monthly_estimates$site <- factor(monthly_estimates$site, 
#                                  levels=c('Redcliffe', 'Toowoomba', 'Sunnybank', 
#                                           'Canungra', 'Burleigh', 'Clunes',
#                                           'Lismore','Nambucca Heads'))
# 
write.csv(monthly_estimates, 'output/monthly_estimates_fieldcases.csv', row.names = FALSE)

#########
hev_prev_loc
hev_prev_loc$date = as.Date(hev_prev_loc$date, format = "%Y-%m-%d")
hev_prev_loc$site_north <- as.factor(hev_prev_loc$northing) #ordered north to south

hev_prev_loc$site = hev_prev_loc$loc
roost_col = wes_palette("Zissou1",8, 'continuous')

# plot comparing monthly estimates to observations
range(hendra_crs$date)
ggplot(monthly_estimates,aes(date, prevalence_est))+
  geom_line(aes(col = site, group = site),  position=position_dodge(0.1))+
  geom_point(data = hev_prev_loc, aes(x = date, y = prev, col = site),  pch=1,position=position_dodge(0.1))+
  # geom_errorbar(data = preempt_shed,aes(x = date_start, y = PointEst, col = site, ymin=Lower, ymax=Upper), width=.2,
  #               position=position_dodge(0.1),
  #               col = 'grey90')+
  theme_test() + 
  scale_y_continuous(limits = c(0,0.8), breaks = seq(0,0.8,0.2))+
  scale_x_date(date_minor_breaks = "1 day", breaks = '4 months', limits = as.Date(c("2011-06-01","2014-12-31")), 
               date_labels = '%Y-%m')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~site, ncol = 3)+
  theme(legend.position = "none")+
  #scale_color_manual(values = roost_col)+
  labs(x = 'date of collection', y = 'prevalence of hendra', col = 'roost')


#######################
# looking at September 2017 observations
trun = preempt_shed[preempt_shed$date_start < as.Date('2017-09-08'),]
trun = trun[trun$date_start > as.Date('2017-08-12'),]
trun
trun_sf = st_as_sf(trun, coords = c( 'easting','northing'),
                   crs = "+proj=utm +zone=56 +south +datum=WGS84")
trun_crs = st_transform(trun_sf, crs_pred)
#plot(prev_field[[10]])
#plot(trun_crs, add = TRUE)
values_buff = terra::extract(prev_field[[9]], buff)$mean

prev_sept = prev_field[[9]]
prev_sept_df = as.data.frame(prev_sept, xy = T)
range(prev_sept_df$mean)

range(trun_crs$PointEst)
ggplot() + 
  geom_sf(data = aus_crs, lwd = 0.5, col = 'black',fill = 'grey90')+
  geom_tile(data = prev_sept_df, aes(x = x,
                                     y = y,
                                     fill = mean), cex = 0.1) +
  scale_fill_viridis_c(limits = c(0, 0.35),
                       n.breaks = 8) +
  labs(x = NULL,
       y = NULL,
       fill = "Expected Prevalence",
       color = "Observed Prevalence") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",
                                    color = "black"),
        panel.grid.major = element_line(color = "grey90")) +
  guides(fill = guide_colorbar(ticks = T,
                               ticks.colour = "black",
                               ticks.linewidth = 1,
                               frame.colour = "black",
                               frame.linewidth = 1,
                               barwidth = 1,
                               barheight = 10))+
  
  geom_sf(data = aus_crs, lwd = 0.5, col = 'black',fill = 'NA')+
  geom_sf(data = trun_crs, aes(col = PointEst), size = 5) +
  geom_sf(data = trun_crs, shape = 1,size = 5,lwd = 2, colour = "white")+
  scale_color_viridis(limits = c(0, 0.35),n.breaks = 8) +
  guides(col = guide_colorbar(ticks = T,
                              ticks.colour = "black",
                              ticks.linewidth = 1,
                              frame.colour = "black",
                              frame.linewidth = 1,
                              barwidth = 1,
                              barheight = 10))+
  
  coord_sf(xlim = c(2636633, 2823360), ylim = c(1996948, 2377938))+
  theme_bw()+
  annotation_scale(location = "bl", width_hint = 0.2) 


