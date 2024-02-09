## Calculate total bacteria abundance and biomass, and mean cell size, in 
## epi, meso and bathypelagic waters

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(terra)
library(RNetCDF)
library(raster)

setwd("~/Desktop/Papers/Bacteria_Census")

bac_abundance <- data.table::fread("contemporary_bacteria.csv")

bathy_data <- read.csv("./stat_model/prediction_data.csv")
bathy_data <- bathy_data %>% dplyr::select(Long, Lat, BATHY)

bacabund_depths <- unique(bac_abundance$depth)

lonlat <- expand.grid("Long" = -179.5:179.5, "Lat" = -89.5:89.5)

### PROCESS NOWICKI RESPIRATION DATA
devries <- open.nc("biopump_model_output.nc")

devries_depths <- var.get.nc(devries, "DEPTH")[1,1,]
bacabund_depths <- unique(bac_abundance$depth)

## Open, reorient, then regrid poc and doc respiration, mmolC/m3/yr
tot_poc_resp <- var.get.nc(devries, "tot_poc_resp")
tot_poc_resp <- aperm(as.array(disagg(rast(tot_poc_resp[c(2:91),c(91:180,1:90),]), 2)), c(2,1,3))
tot_doc_resp <- var.get.nc(devries, "tot_doc_resp")
tot_doc_resp <- aperm(as.array(disagg(rast(tot_doc_resp[c(2:91),c(91:180,1:90),]), 2)), c(2,1,3))
tot_all_resp <- var.get.nc(devries, "tot_resp")
tot_all_resp <- aperm(as.array(disagg(rast(tot_all_resp[c(2:91),c(91:180,1:90),]), 2)), c(2,1,3))
tot_zoo_resp <- var.get.nc(devries, "zoo_resp")
tot_zoo_resp <- aperm(as.array(disagg(rast(tot_zoo_resp[c(2:91),c(91:180,1:90),]), 2)), c(2,1,3))
tot_vol <- var.get.nc(devries, "Volume")

#tot_npp <- apply(var.get.nc(devries, "NPP"), c(1,2), mean, na.rm = TRUE)
#tot_npp <- aperm(as.array(disagg(rast(tot_npp[c(2:91),c(91:180,1:90)]), 2)), c(2,1,3))
#tot_npp <- sum(tot_npp[,,1]*t(as.matrix(raster::area(raster()))), na.rm = TRUE)*12.01/(1000*1e15)*1e6

#tot_exp_carb <- var.get.nc(devries, "tot_ex")
#tot_exp_carb <- aperm(as.array(disagg(rast(tot_exp_carb[c(2:91),c(91:180,1:90)]), 2)), c(2,1,3))
#tot_exp_carb <- sum(tot_exp_carb[,,1]*t(as.matrix(raster::area(raster()))), na.rm = TRUE)*12.01/(1000*1e15)*1e6

tot_poc_resp2 <- array(NA, dim = c(360, 180, length(bacabund_depths)))
tot_doc_resp2 <- array(NA, dim = c(360, 180, length(bacabund_depths)))
tot_all_resp2 <- array(NA, dim = c(360, 180, length(bacabund_depths)))
tot_zoo_resp2 <- array(NA, dim = c(360, 180, length(bacabund_depths)))

for(i in 1:length(devries_depths)){
  if(i == 1){
    curr_idx <-which(bacabund_depths <= devries_depths[i])
  }
  
  if(i > 1){
    curr_idx <-which(bacabund_depths <= devries_depths[i] & bacabund_depths > devries_depths[i-1]) 
  }
  
  if(i == length(devries_depths)){
    curr_idx <- c(curr_idx, length(bacabund_depths))
  }
  
  tot_poc_resp2[,,curr_idx] <- tot_poc_resp[,,i]
  tot_doc_resp2[,,curr_idx] <- tot_doc_resp[,,i]
  tot_all_resp2[,,curr_idx] <- tot_all_resp[,,i]
  tot_zoo_resp2[,,curr_idx] <- tot_zoo_resp[,,i]
}

lonlat <- expand.grid("Long" = -179.5:179.5, "Lat" = -89.5:89.5)

resp_data <- data.frame("Long" = rep(lonlat$Long, times = length(bacabund_depths)), "Lat" = rep(lonlat$Lat, times = length(bacabund_depths)), 
                        "Depth" = rep(bacabund_depths, each = dim(lonlat)[1]),"POC_resp_mmolC_m3_yr" = as.vector(tot_poc_resp2), "DOC_resp_mmolC_m3_yr" = as.vector(tot_doc_resp2),
                        "All_resp_mmolC_m3_yr" = as.vector(tot_all_resp2),
                        "Zoo_resp_mmolC_m3_yr" = as.vector(tot_zoo_resp2))

resp_data <- resp_data %>% dplyr::mutate(
  Delta_resp = All_resp_mmolC_m3_yr - Zoo_resp_mmolC_m3_yr - POC_resp_mmolC_m3_yr - DOC_resp_mmolC_m3_yr)

kk <- resp_data %>% dplyr::filter(Depth == 0.5)

bac_abundance <- left_join(bac_abundance, resp_data, by = c("Long" = "Long", "Lat" = "Lat", "depth" = "Depth"))

bac_abundance <- bac_abundance %>% dplyr::select(Long, Lat, logdepth, log_nitrate, temperature, aou, logchlo, 
                                                 log_bacabund_ml, bac_cell_carbon_fgC, bac_biomass_gCm3, POC_resp_mmolC_m3_yr, 
                                                 DOC_resp_mmolC_m3_yr, All_resp_mmolC_m3_yr, Zoo_resp_mmolC_m3_yr)
bac_abundance <- bac_abundance %>% rename(log10_depth_m = logdepth,
                                          log10_nitrate_umol_kg = log_nitrate,
                                          log10_chlo_mg_m3 = logchlo,
                                          aou_umol_kg = aou,
                                          temperature_C = temperature,
                                          log10_bacabund_per_ml = log_bacabund_ml) %>%
                                    dplyr::mutate(bac_resp_gC_m3_d = 12.01/(1000*365)*(POC_resp_mmolC_m3_yr + DOC_resp_mmolC_m3_yr),
                                                  bac_resp_d = bac_resp_gC_m3_d*1e15/((10^log10_bacabund_per_ml)*1e6*bac_cell_carbon_fgC),
                                                  all_resp_gC_m3_d = 12.01/(1000*365)*All_resp_mmolC_m3_yr,
                                                  zoo_resp_gC_m3_d = 12.01/(1000*365)*Zoo_resp_mmolC_m3_yr) %>%
                        left_join(bathy_data, by = c("Long" = "Long", "Lat" = "Lat")) %>% ungroup()


npp_data <- read.csv("./npp_data/intpp_clim.csv")

npp_mean <- npp_data  %>% rowwise() %>% 
  mutate(intpp_gC_m2_d=mean(c(Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec), na.rm = TRUE)/1000) %>% dplyr::select(lon, lat, intpp_gC_m2_d)
npp_mean$area_m2 <- 1e6*as.vector(as.matrix(area(raster())))

## Extract top 200m (euphotic zone)
euph <- bac_abundance %>% filter(log10_depth_m <= log10(200)) %>% group_by(Long, Lat) %>%
  dplyr::mutate(sgr_d = 10^(-1.766477+0.0311771*temperature_C+0.134835*log10_chlo_mg_m3), # SGR, from data
                total_growth_gC_m3_d = sgr_d*bac_biomass_gCm3,
                total_resp_gC_m3_d = total_growth_gC_m3_d*(1-0.14)/0.14, # Convert growth to respiration, assuming a constant GGE of 0.14 (Pep data)
                total_resp_urrutita_gC_m3_d = 0.1*(10^log10_bacabund_per_ml*1e6*1e-15*3.21e11*exp(-0.589/(8.62e-5*(273.15+temperature_C)))), # Lopez-Urrutia 2007 equation, assuming 10% of prokaryotes are active
                depth_m = 10^(log10_depth_m),
                depth_m = if_else(depth_m < max(depth_m), depth_m, 200), # Stretch deepest interval to go to 200m
                depth_m = ifelse(depth_m > BATHY, BATHY, depth_m), # If deepest interval bottom exceeds bathy, swap to bathy
                depth_abs = depth_m - dplyr::lag(depth_m, default=0)) %>%
  dplyr::summarise(temperature_C = sum(temperature_C*depth_abs, na.rm = TRUE)/sum(depth_abs, na.rm = TRUE), # mean temperature 
                   aou_umol_kg = sum(aou_umol_kg*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean aou
                   nitrate_umol_kg = sum(10^log10_nitrate_umol_kg*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean nitrate
                   bacabund_m2 = sum(10^log10_bacabund_per_ml*1e6*depth_abs, na.rm = TRUE), # total bacteria abundance
                   bacabund_m3 = sum(10^log10_bacabund_per_ml*1e6*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean bacteria abundance
                   fgC_cell = sum(bac_cell_carbon_fgC*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean bacteria cell size
                   bac_biomass_gC_m2 = sum(bac_biomass_gCm3*depth_abs, na.rm = TRUE),# total bacteria carbon biomass
                   all_resp_gC_m2_d = sum(all_resp_gC_m3_d*depth_abs, na.rm = TRUE), # all respiration
                   bac_biomass_gC_m3 = sum(bac_biomass_gCm3*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean bacteria biomass
                   bac_resp_gC_m2_d = sum(bac_resp_gC_m3_d*depth_abs, na.rm = TRUE), # total bacteria respiration
                   bac_resp_d = bac_resp_gC_m2_d/bac_biomass_gC_m2,
                   all_zoo_gC_m2_d = sum(zoo_resp_gC_m3_d*depth_abs, na.rm = TRUE), # all respiration
                   mean_sgr_d = sum(sgr_d*depth_abs, na.rm = TRUE)/sum(depth_abs, na.rm = TRUE),
                   bac_resp2_gC_m2_d = sum(total_resp_gC_m3_d*depth_abs, na.rm = TRUE), # total bacteria respiration (from data)
                   total_resp_urrutia_gC_m2_d = sum(total_resp_urrutita_gC_m3_d*depth_abs, na.rm = TRUE),
                   bac_growth2_gC_m2_d = sum(total_growth_gC_m3_d*depth_abs, na.rm = TRUE), # total bacteria growth (from data) 
                   bac_c_demand_gC_m2_d = bac_resp2_gC_m2_d + bac_growth2_gC_m2_d, # total bacteria carbon demand (from data)
                   bac_c_demand_gC_m3_d = (bac_resp2_gC_m2_d + bac_growth2_gC_m2_d)/sum(depth_abs), # mean total bacteria carbon demand
                   .groups = "drop") %>% # specific bacteria respiration
                  left_join(npp_mean, by = c("Long" = "lon", "Lat" = "lat"))

ggplot(euph) + geom_raster(aes(x=Long,y=Lat,fill=bac_resp2_gC_m2_d/intpp_gC_m2_d)) + scale_fill_gradientn(limits = c(0,1), colours = rainbow(12))



### Extract top 200m (epipelagic)
epip <- bac_abundance %>% filter(log10_depth_m <= log10(200)) %>% group_by(Long, Lat) %>%
  dplyr::mutate(sgr_d = 10^(-1.766477+0.0311771*temperature_C+0.134835*log10_chlo_mg_m3), # SGR, from data
                total_growth_gC_m3_d = sgr_d*bac_biomass_gCm3,
                gge1 = 1/((0.037 + 0.65*(total_growth_gC_m3_d*1e3/24))/(1.8 + total_growth_gC_m3_d*1e3/24)), # from del Giorgio and Cole, 1998
                gge2 = 1/(1-1/(0.727*10^log10_chlo_mg_m3/(10^log10_chlo_mg_m3 + 4.08) + 1.02)), # from Lopez-Urrutia 2007
                total_resp_gC_m3_d = total_growth_gC_m3_d*(1-0.14)/0.14, # Convert growth to respiration, assuming a constant GGE of 0.14 (from Pep data)
                depth_m = 10^(log10_depth_m),
                depth_m = if_else(depth_m < max(depth_m), depth_m, 200), # Stretch deepest interval to go to 200m
                depth_m = ifelse(depth_m > BATHY, BATHY, depth_m), # If deepest interval bottom exceeds bathy, swap to bathy
                depth_abs = depth_m - dplyr::lag(depth_m, default=0))%>%
  dplyr::summarise(temperature_C = sum(temperature_C*depth_abs, na.rm = TRUE)/sum(depth_abs, na.rm = TRUE), # mean temperature 
                   aou_umol_kg = sum(aou_umol_kg*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean aou
                   nitrate_umol_kg = sum(10^log10_nitrate_umol_kg*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean nitrate
                   bacabund_m2 = sum(10^log10_bacabund_per_ml*1e6*depth_abs, na.rm = TRUE), # total bacteria abundance
                   bacabund_m3 = sum(10^log10_bacabund_per_ml*1e6*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean bacteria abundance
                   fgC_cell = sum(bac_cell_carbon_fgC*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean bacteria cell size
                   bac_biomass_gC_m2 = sum(bac_biomass_gCm3*depth_abs, na.rm = TRUE),# total bacteria carbon biomass
                   all_resp_gC_m2_d = sum(all_resp_gC_m3_d*depth_abs, na.rm = TRUE), # all respiration
                   bac_biomass_gC_m3 = sum(bac_biomass_gCm3*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean bacteria biomass
                   bac_resp_gC_m2_d = sum(bac_resp_gC_m3_d*depth_abs, na.rm = TRUE), # total bacteria respiration
                   bac_resp_d = bac_resp_gC_m2_d/bac_biomass_gC_m2,
                   all_zoo_gC_m2_d = sum(zoo_resp_gC_m3_d*depth_abs, na.rm = TRUE), # all respiration
                   mean_sgr_d = sum(sgr_d*depth_abs, na.rm = TRUE)/sum(depth_abs, na.rm = TRUE),
                   bac_resp2_gC_m2_d = sum(total_resp_gC_m3_d*depth_abs, na.rm = TRUE), # total bacteria respiration (from data)
                   bac_growth2_gC_m2_d = sum(total_growth_gC_m3_d*depth_abs, na.rm = TRUE), # total bacteria growth (from data) 
                   bac_c_demand_gC_m2_d = bac_resp2_gC_m2_d + bac_growth2_gC_m2_d, # total bacteria carbon demand (from data)
                   bac_c_demand_gC_m3_d = (bac_resp2_gC_m2_d + bac_growth2_gC_m2_d)/sum(depth_abs), # mean total bacteria carbon demand
                   .groups = "drop") %>% # specific bacteria respiration
                  left_join(npp_mean, by = c("Long" = "lon", "Lat" = "lat"))

epip <- left_join(lonlat, epip, by = c("Long" = "Long", "Lat" = "Lat")) # Add in NA cells

ggplot(epip) + geom_raster(aes(x=Long,y=Lat,fill=bac_resp2_gC_m2_d/intpp_gC_m2_d)) + scale_fill_gradientn(limits = c(0,5), colours = rainbow(12))

epip$area_m2 <- 1e6*as.vector(t(as.matrix(area(raster()))))
sum(epip$area_m2*epip$bac_resp2_gC_m2_d*365, na.rm = TRUE)/1e15

write.csv(epip, "bacteria_epipelagic.csv", row.names = FALSE)


ggplot(epip) + geom_raster(aes(x=Long,y=Lat,fill=bac_resp2_gC_m2_d/intpp_gC_m2_d)) + 
  scale_fill_gradientn(limits = c(0,1), colours = rainbow(12), name = "Ratio") +
  labs(subtitle = "Prok. Respiration/Intpp")

ggplot(epip) + geom_raster(aes(x=Long,y=Lat,fill=bac_c_demand_gC_m2_d/intpp_gC_m2_d)) + 
  scale_fill_gradientn(limits = c(0,1), colours = rainbow(12), name = "Ratio") +
  labs(subtitle = "Prok. Total/Intpp")

### Extract 200-1000m (mesopelagic)
mesop <- bac_abundance %>% dplyr::filter(log10_depth_m <= log10(1000) & log10_depth_m > log10(200))%>% group_by(Long, Lat) %>% 
  dplyr::mutate(depth_m = 10^(log10_depth_m),
                depth_m = ifelse(depth_m == min(depth_m), 200, depth_m), # Stretch top of shallowest interval to 200m
                depth_m = if_else(depth_m < max(depth_m), depth_m, 1000), # Stretch deepest interval to go to 200m
                depth_m = ifelse(depth_m > BATHY, BATHY, depth_m), # If deepest interval bottom exceeds bathy, swap to bathy
                depth_abs = depth_m - dplyr::lag(depth_m, default=0))%>%
  dplyr::summarise(temperature_C = sum(temperature_C*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean temperature 
                   aou_umol_kg = sum(aou_umol_kg*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean aou
                   nitrate_umol_kg = sum(10^log10_nitrate_umol_kg*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean nitrate
                   bacabund_m2 = sum(10^log10_bacabund_per_ml*1e6*depth_abs, na.rm = TRUE), # total bacteria abundance
                   bacabund_m3 = sum(10^log10_bacabund_per_ml*1e6*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean bacteria abundance
                   fgC_cell = sum(bac_cell_carbon_fgC*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean bacteria cell size
                   bac_biomass_gC_m2 = sum(bac_biomass_gCm3*depth_abs, na.rm = TRUE),# total bacteria carbon biomass
                   bac_biomass_gC_m3 = sum(bac_biomass_gCm3*depth_abs, na.rm = TRUE)/sum(depth_abs),
                   bac_resp_gC_m2_d = sum(bac_resp_gC_m3_d*depth_abs, na.rm = TRUE), # total bacteria respiration
                   bac_resp_d = bac_resp_gC_m2_d/bac_biomass_gC_m2,
                   all_resp_gC_m2_d = sum(all_resp_gC_m3_d*depth_abs, na.rm = TRUE), # all respiration
                   all_zoo_gC_m2_d = sum(zoo_resp_gC_m3_d*depth_abs, na.rm = TRUE), 
                   .groups = "drop") # all respiration 

mesop <- left_join(lonlat, mesop, by = c("Long" = "Long", "Lat" = "Lat"))
write.csv(mesop, "bacteria_mesopelagic.csv", row.names = FALSE)
  
### Extract >1000m (bathyelagic)
bathyp <- bac_abundance %>% dplyr::filter(log10_depth_m >= log10(1000))%>% group_by(Long, Lat) %>% 
  dplyr::mutate(depth_m = 10^(log10_depth_m),
                depth_m = ifelse(depth_m == min(depth_m), 1000, depth_m), # Stretch top of shallowest interval to 200m
                depth_m = if_else(depth_m < BATHY, depth_m, BATHY), # Stretch deepest interval to go to sea floor
                depth_m = ifelse(depth_m > BATHY, BATHY, depth_m), # If deepest interval bottom exceeds bathy, swap to bathy
                depth_abs = depth_m - dplyr::lag(depth_m, default=0))%>%
  dplyr::summarise(temperature_C = sum(temperature_C*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean temperature 
                   aou_umol_kg = sum(aou_umol_kg*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean aou
                   nitrate_umol_kg = sum(10^log10_nitrate_umol_kg*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean nitrate
                   bacabund_m2 = sum(10^log10_bacabund_per_ml*1e6*depth_abs, na.rm = TRUE), # total bacteria abundance
                   bacabund_m3 = sum(10^log10_bacabund_per_ml*1e6*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean bacteria abundance
                   fgC_cell = sum(bac_cell_carbon_fgC*depth_abs, na.rm = TRUE)/sum(depth_abs), # mean bacteria cell size
                   bac_biomass_gC_m2 = sum(bac_biomass_gCm3*depth_abs, na.rm = TRUE),# total bacteria carbon biomass
                   bac_biomass_gC_m3 = sum(bac_biomass_gCm3*depth_abs, na.rm = TRUE)/sum(depth_abs),
                   bac_resp_gC_m2_d = sum(bac_resp_gC_m3_d*depth_abs, na.rm = TRUE), # total bacteria respiration
                   bac_resp_d = bac_resp_gC_m2_d/bac_biomass_gC_m2, # specific bacteria respiration 
                   all_resp_gC_m2_d = sum(all_resp_gC_m3_d*depth_abs, na.rm = TRUE), # all respiration
                   all_zoo_gC_m2_d = sum(zoo_resp_gC_m3_d*depth_abs, na.rm = TRUE), 
                   .groups = "drop") # all respiration 

bathyp <- left_join(lonlat, bathyp, by = c("Long" = "Long", "Lat" = "Lat"))
write.csv(bathyp, "bacteria_bathypelagic.csv", row.names = FALSE)

ggplot(epip) + geom_raster(aes(x=Long,y=Lat,fill=temperature_C))
ggplot(mesop) + geom_raster(aes(x=Long,y=Lat,fill=temperature_C))
ggplot(bathyp) + geom_raster(aes(x=Long,y=Lat,fill=temperature_C))
