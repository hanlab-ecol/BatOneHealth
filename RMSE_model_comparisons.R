# RMSE Model comparisons

# libraries
library(terra)

#model output downloaded from here: https://drive.google.com/drive/folders/15OiPaN1dYh3NIBJU2ne3f0srw7uZ8e8o
loc_dir = '/Users/cfaust/Documents/workspace/hendra_msm_large_files/Results_20220811'
rast = lapply(list.files(path = file.path(loc_dir,'0 Averaging/Prevalence'), pattern = "^E", full.names = T), terra::rast)
