# RMSE Model comparisons

# libraries
library(terra)
library(dplyr)
library(sf)
library(wesanderson)
library(ggplot2)
library(viridis)
library(ggspatial)
library(zoo)
library(hydroGOF)

################
# importing Huma Fields data
hev_raw = read.csv('data/hendra-virus-transmission-data-east-australia.csv', header = TRUE)
hev_loc = read.csv('data/locations_sites.csv', header = TRUE)
range(hev_prev_loc$UTM_zone)
range(hev_prev_loc$northing)
head(hev_raw)
hev_raw$date = as.Date(hev_raw$date, format = '%d/%m/%Y')
hev_prev = as.data.frame(hev_raw %>%
                           group_by(date,loc) %>%
                           summarize(prev = mean(res, na.rm = TRUE),
                                     obs = length(res),
                                     fftotal = mean(fftotal, na.rm = TRUE),
                                     estb = mean(estb, na.rm = TRUE),
                                     pct_b = mean(pct_b, na.rm = TRUE)))

hev_prev_MEAN = as.data.frame(hev_raw %>%
                                group_by(loc) %>%
                                summarize(prev = mean(res, na.rm = TRUE),
                                          obs = length(res),
                                          fftotal = mean(fftotal, na.rm = TRUE),
                                          estb = mean(estb, na.rm = TRUE),
                                          pct_b = mean(pct_b, na.rm = TRUE)))

hev_prev_loc = merge(hev_prev, hev_loc, by.x = 'loc', by.y ='site')
hev_prev_MEAN_loc = merge(hev_prev_MEAN, hev_loc, by.x = 'loc', by.y ='site') #per site...

#write.csv(hev_prev, "output/hev_prev_summary_field.csv", row.names = FALSE)
#write.csv(hev_prev_MEAN_loc, "output/hev_prev_summary_bylocation_wloc_field.csv", row.names = FALSE)

ggplot(hev_prev, aes(x=prev,pct_b))+
  geom_point()


##############
# partitioning into different UTMs and turing into a spatial DF
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
# australia shapefile
ausdir = "~/Dropbox/Projects/Hendra/workspace/sp_dat/AUS_adm/nsaasr9nnd_02211a04es_geo___"
aus = st_read(ausdir,'aust_cd66states', crs = "+proj=longlat +datum=WGS84") #country boundary
aus_crs = st_transform(aus, crs_pred)
aus_simp <- st_simplify(aus_crs, preserveTopology = FALSE, 
                        dTolerance = 10000)


buff = st_buffer(hendra_crs[1,], dist = 5000)
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

range(hendra_crs$date)

###############
# importing prevalence rasters and transforming
all_dir = "/Users/cfaust/Documents/workspace/hendra_msm_large_files/Results_20220811/Models/AllModels"
all_temp <- list.files(path = all_dir,pattern = "^E", full.names = TRUE)
all_rast <- rast(all_temp[5:8]) # years overlapping with shedding data

foodshort_dir = "/Users/cfaust/Documents/workspace/hendra_msm_large_files/Results_20220811/Models/FoodShortage"
foodshort_temp <- list.files(path = foodshort_dir,pattern = "^E", full.names = TRUE)
foodshort_rast <- rast(foodshort_temp[5:8]) # years overlapping with shedding data

foodshortnew_dir = "/Users/cfaust/Documents/workspace/hendra_msm_large_files/Results_20220811/Models/FoodShortageNewRoost"
foodshortnew_temp <- list.files(path = foodshortnew_dir,pattern = "^E", full.names = TRUE)
foodshortnew_rast <- rast(foodshortnew_temp[5:8]) # years overlapping with shedding data

foodshortrehab_dir = "/Users/cfaust/Documents/workspace/hendra_msm_large_files/Results_20220811/Models/FoodShortageRehab"
foodshortrehab_temp <- list.files(path = foodshortrehab_dir,pattern = "^E", full.names = TRUE)
foodshortrehab_rast <- rast(foodshortrehab_temp[5:8]) # years overlapping with shedding data

new_dir = "/Users/cfaust/Documents/workspace/hendra_msm_large_files/Results_20220811/Models/NewRoost"
new_temp <- list.files(path = new_dir,pattern = "^E", full.names = TRUE)
new_rast <- rast(new_temp[5:8]) # years overlapping with shedding data

newrehab_dir = "/Users/cfaust/Documents/workspace/hendra_msm_large_files/Results_20220811/Models/NewRoostRehab"
newrehab_temp <- list.files(path = newrehab_dir,pattern = "^E", full.names = TRUE)
newrehab_rast <- rast(newrehab_temp[5:8]) # years overlapping with shedding data

rehab_dir = "/Users/cfaust/Documents/workspace/hendra_msm_large_files/Results_20220811/Models/Rehab"
rehab_temp <- list.files(path = rehab_dir,pattern = "^E", full.names = TRUE)
rehab_rast <- rast(rehab_temp[5:8]) # years overlapping with shedding data





monthly_estimates = as.data.frame(matrix(NA, nrow = nrow(hendra_crs)*nlyr(all_rast), ncol =7))
colnames(monthly_estimates) = c('all_exp', 'foodshortage_exp','foodshortagenew_exp','foodshortagerehab_exp',
                                'new_exp','newrehab_exp','rehab_exp')
monthly_estimates$site = rep(hendra_crs$loc, each=nlyr(all_rast))

for(j in 301:nrow(hendra_crs)){#nrow(hendra_crs)
  #j = 6
  buff =  vect(st_buffer(hendra_crs[j,], dist = 1000))
  for(i in 1:nlyr(all_rast)){
    #i = 10
    #ALL
    values_buff = terra::extract(all_rast[[i]], buff)$mean
    monthly_estimates$all_exp[(j-1)*nlyr(all_rast)+i] = mean(values_buff[!is.na(values_buff)]) 
    #food short
    values_buff2 = terra::extract(foodshort_rast[[i]], buff)$mean
    monthly_estimates$foodshortage_exp[(j-1)*nlyr(all_rast)+i] = mean(values_buff2[!is.na(values_buff2)]) 
    #food short new
    values_buff3 = terra::extract(foodshortnew_rast[[i]], buff)$mean
    monthly_estimates$foodshortagenew_exp[(j-1)*nlyr(all_rast)+i] = mean(values_buff3[!is.na(values_buff3)]) 
    #food short rehab
    values_buff4 = terra::extract(foodshortrehab_rast[[i]], buff)$mean
    monthly_estimates$foodshortagerehab_exp[(j-1)*nlyr(all_rast)+i] = mean(values_buff4[!is.na(values_buff4)]) 
    
    values_buff5 = terra::extract(new_rast[[i]], buff)$mean
    monthly_estimates$new_exp[(j-1)*nlyr(all_rast)+i] = mean(values_buff5[!is.na(values_buff5)]) 
    
    values_buff6 = terra::extract(newrehab_rast[[i]], buff)$mean
    monthly_estimates$newrehab_exp[(j-1)*nlyr(all_rast)+i] = mean(values_buff6[!is.na(values_buff6)]) 
    
    values_buff7 = terra::extract(rehab_rast[[i]], buff)$mean
    monthly_estimates$rehab_exp[(j-1)*nlyr(all_rast)+i] = mean(values_buff7[!is.na(values_buff7)]) 
    
  }
}
head(monthly_estimates)
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
write.csv(monthly_estimates, 'output/monthly_estimates_fieldcases_20220818.csv', row.names = FALSE)

#########
hev_prev_loc
hev_prev_loc$date = as.Date(hev_prev_loc$date, format = "%Y-%m-%d")
hev_prev_loc$site_north <- as.factor(hev_prev_loc$northing) #ordered north to south

hev_prev_loc$site = hev_prev_loc$loc
roost_col = wes_palette("Zissou1",8, 'continuous')

# plot comparing monthly estimates to observations
range(hendra_crs$date)
ggplot(monthly_estimates,aes(date, all_exp))+
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

###########
monthly_estimates$monyr = as.yearmon(monthly_estimates$date)
hev_prev_loc$monyr = as.yearmon(hev_prev_loc$date)

monthly_est_obs = inner_join(monthly_estimates, hev_prev_loc, by = c('monyr', 'site'))
monthly_est_obs2 = monthly_est_obs[!duplicated(monthly_est_obs),]
rmse = sqrt(mean((monthly_est_obs2$prev - monthly_est_obs2$all_exp)^2, na.rm = TRUE))


#It looks like we should be normalizing RMSE by either dividing by the SD of the observed values or the Max - Min of the observed values. I am seeing it pop up both ways and haven't found a foundational explanation of one vs. the other yet. One R package, 'hydroGOF', calculates NRMSE with SD as the default and max-min as the alternative option
head(monthly_est_obs2)
rmse_roost = as.data.frame(monthly_est_obs2 %>%
                             group_by(site) %>%
                             summarize(rmse_all = sqrt(mean((prev - all_exp)^2, na.rm=TRUE)),
                                       rmse_foodshortage = sqrt(mean((prev - foodshortage_exp)^2, na.rm=TRUE)),
                                       rmse_foodshortagenew = sqrt(mean((prev - foodshortagenew_exp)^2, na.rm=TRUE)),
                                       rmse_foodshortagerehab = sqrt(mean((prev - foodshortagerehab_exp)^2, na.rm=TRUE)),
                                       rmse_new = sqrt(mean((prev - new_exp)^2, na.rm=TRUE)),
                                       rmse_newrehab = sqrt(mean((prev - newrehab_exp)^2, na.rm=TRUE)),
                                       rmse_rehab = sqrt(mean((prev - rehab_exp)^2, na.rm=TRUE)),
                                       
                                       nrmse_all = nrmse(all_exp, prev, na.rm=TRUE, norm="sd"),
                                       nrmse_foodshortage = nrmse(foodshortage_exp, prev, na.rm=TRUE, norm="sd"),
                                       nrmse_foodshortagenew = nrmse(foodshortagenew_exp, prev, na.rm=TRUE, norm="sd"),
                                       nrmse_foodshortagerehab = nrmse(foodshortagerehab_exp, prev, na.rm=TRUE, norm="sd"),
                                       nrmse_new = nrmse(new_exp, prev, na.rm=TRUE, norm="sd"),
                                       nrmse_newrehab = nrmse(newrehab_exp, prev, na.rm=TRUE, norm="sd") ,
                                       nrmse_rehab = nrmse(rehab_exp, prev, na.rm=TRUE, norm="sd"),
                                       
                                       nrmse2_all = nrmse(all_exp, prev, na.rm=TRUE, norm="maxmin"),
                                       nrmse2_foodshortage = nrmse(foodshortage_exp, prev, na.rm=TRUE, norm="maxmin"),
                                       nrmse2_foodshortagenew = nrmse(foodshortagenew_exp, prev, na.rm=TRUE, norm="maxmin"),
                                       nrmse2_foodshortagerehab = nrmse(foodshortagerehab_exp, prev, na.rm=TRUE, norm="maxmin"),
                                       nrmse2_new = nrmse(new_exp, prev, na.rm=TRUE, norm="maxmin"),
                                       nrmse2_newrehab = nrmse(newrehab_exp, prev, na.rm=TRUE, norm="maxmin") ,
                                       nrmse2_rehab = nrmse(rehab_exp, prev, na.rm=TRUE, norm="maxmin"),
                                       
                                       
                                       obs = length(prev),
                                       mean_obs = mean(prev),
                                       mean_all = mean(all_exp, na.rm = TRUE),
                                       mean_foodshortage= mean(foodshortage_exp, na.rm = TRUE),
                                       mean_foodshortagenew= mean(foodshortagenew_exp, na.rm = TRUE),
                                       mean_foodshortagerehab= mean(foodshortagerehab_exp, na.rm = TRUE),
                                       mean_new= mean(new_exp, na.rm = TRUE),
                                       mean_newrehab= mean(newrehab_exp, na.rm = TRUE),
                                       mean_rehab= mean(rehab_exp, na.rm = TRUE)
                             ))

rmse_roost
#normalized between each roost
#overall measure > roosts with more data points
cor(rmse_roost$mean_all, rmse_roost$mean_newrehab) #0.9844378
cor(rmse_roost$mean_all, rmse_roost$mean_foodshortage) # -0.2757267
plot(rmse_roost$mean_all, rmse_roost$mean_foodshortage) # -0.2757267

cor(rmse_roost$mean_all, rmse_roost$mean_foodshortagenew) # 0.2944876
plot(rmse_roost$mean_all, rmse_roost$mean_foodshortagenew)

cor(rmse_roost$mean_all, rmse_roost$mean_foodshortagerehab) # 0.9717971
plot(rmse_roost$mean_all, rmse_roost$mean_foodshortagerehab)

cor(rmse_roost$mean_all, rmse_roost$mean_new) # 0.2268154
plot(rmse_roost$mean_all, rmse_roost$mean_new) 

cor(rmse_roost$mean_all, rmse_roost$mean_rehab) # 0.9739941
plot(rmse_roost$mean_all, rmse_roost$mean_rehab)


#new rehab

plot(monthly_est_obs2$prev, monthly_est_obs2$all_exp)



