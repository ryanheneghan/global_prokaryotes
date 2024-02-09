rm(list = ls())
library(RNetCDF)
library(dplyr)
library(tidyr)
library(data.table)
library(raster)

setwd("~/Desktop/Papers/Bacteria_Census")

## For fish biomass delta, use Lotze et al. 2019 relationship of 5% decline in global
## fish biomass per 1C warming.

## For fish respiration delta, use Lopez-Urrutia et al. 2006 activation energy of -0.56 for
## fish (heterotroph) metabolism

load("microzoo_model_hatton.RData") # microzoo biomass model
microzoo_model <- gm1
load("rhizaria_top200_model_hatton.RData") # rhizaria biomass model
rhizaria_model <- gm1
load("macrozoo_model_hatton.RData") # macrozoo biomass model
macrozoo_model <- gm1
load("mesozoo_model_hatton.RData") # mesozoo biomass model
mesozoo_model <- gm1

load("lm_model_final.RData") # bacteria abundance model
load("glmm_model1_carbon.RData") # bacteria cell carbon model


bac_contemp <- data.table::fread("bacteria_epipelagic.csv")
bac_contemp <- bac_contemp %>% dplyr::filter(bac_resp_d != 0) %>% dplyr::rename(temperature = temperature_C)

fish_contemp <- data.table::fread("./biomass_other_heterotrophs/fish_predictions.csv")

bathy <- read.csv("./npp_data/prediction_data.csv")
bathy <- bathy %>% dplyr::select(Long, Lat, BATHY)
area_grid <- t(as.matrix(raster::area(raster())))
bathy$area_m2 <- as.vector(area_grid)*1e6*bathy$BATHY/bathy$BATHY

fish_contemp <- left_join(fish_contemp, bathy, by = c("lon" = "Long", "lat" = "Lat"))
bac_contemp <- left_join(bac_contemp, bathy, by = c("Long" = "Long", "Lat" = "Lat"))

bathy$area_km2 <- as.vector(area_grid)

esms <- c("GFDL-ESM4", "IPSL-CM6A-LR", "CMCC-ESM2", "MPI-ESM1-2-LR")
clims <- c("historical", "ssp370")
start_year <- c(1950, 2015, 2015, 2015)
end_year <- c(2014, 2100, 2100, 2100)
depths <- c("lev", "olevel","lev", "lev")
esm_scalars <- read.csv("./esm_forcings/esm_scalars.csv")

### STEP 1
## Save output from each esm, one at a time (memory limit reached otherwise)

for(i in 1:length(esms)){
  curr_esm <- esms[i] 
  print(curr_esm)
  
  for(j in 1:length(clims)){
    curr_clim <- clims[j]
    print(curr_clim)
    ## THETAO
    curr_clim_file <- open.nc(list.files(path = paste("./esm_forcings/", curr_esm, "/", sep = ""), pattern = glob2rx(paste("thetao*", curr_clim, "*",sep = "")), full.names = TRUE))
    theta_correct <- esm_scalars[2,i+1]
    
    curr_depths <- var.get.nc(curr_clim_file, depths[i])
    
    depth_idx <- which(curr_depths < 200)
    
    if(length(depth_idx) > 4){ # IPSL AND CMCC HAVE 11 AND 8 DEPTH LEVELS <200M, TOO BIG TO PROCESS, SO TRIM
      depth_idx <- depth_idx[seq(1,length(depth_idx),3)]
    }
    
    curr_depths <- curr_depths[depth_idx]
    
    curr_thetao_file <- var.get.nc(curr_clim_file, "thetao")
    curr_thetao <- curr_thetao_file[c(181:360,1:180),,depth_idx,] + theta_correct
    rm(curr_thetao_file)
    
    ## SST
    SST <- curr_thetao[,,1,] # Extract surface values
    SST <- aperm(replicate(dim(curr_thetao)[3],SST), c(1,2,4,3)) # Overwrite all non-surface values with surface values (lazy)
    
    
    ## PHYC
    curr_phyc_file <- open.nc(list.files(path = paste("./esm_forcings/", curr_esm, "/", sep = ""), pattern = glob2rx(paste("phyc*", curr_clim, "*",sep = "")), full.names = TRUE))

    curr_depths <- var.get.nc(curr_phyc_file, depths[i])
    
    depth_idx <- which(curr_depths < 200)
    
    if(length(depth_idx) > 4){ # IPSL AND CMCC HAVE 11 AND 8 DEPTH LEVELS <200M, TOO BIG TO PROCESS, SO TRIM
      depth_idx <- depth_idx[seq(1,length(depth_idx),3)]
    }
    
    curr_depths <- curr_depths[depth_idx]
    depth_integrals <- diff(c(0,curr_depths[-1], 200)) # Depth intervals in top 200m
    
    curr_phyc_file <- var.get.nc(curr_phyc_file, "phyc")
    phyc <- curr_phyc_file[c(181:360,1:180),,depth_idx,]

    rm(curr_phyc_file)
    
    ## INTPP
    curr_intpp_file <- open.nc(list.files(path = paste("./esm_forcings/", curr_esm, "/", sep = ""), pattern = glob2rx(paste("intpp*", curr_clim, "*",sep = "")), full.names = TRUE))
    
    intpp <- var.get.nc(curr_intpp_file, "intpp")
    intpp <- aperm(replicate(length(depth_idx), intpp[c(181:360,1:180),,]), c(1,2,4,3)) # Overwrite all non-surface values with surface values (lazy)

    rm(curr_intpp_file)
    
    ## CHL
    curr_chl_file <- var.get.nc(open.nc(list.files(path = paste("./esm_forcings/", curr_esm, "/", sep = ""), pattern = glob2rx(paste("chl*", curr_clim, "*",sep = "")), full.names = TRUE)), "chl")
    chl_correct <- esm_scalars[1,i+1]
    
    if(curr_esm == "IPSL-CM6A-LR"){
      curr_chl <- curr_chl_file[c(181:360,1:180),,,]*1e3*chl_correct ## Extract surface, convert from g m-3 to mg m-3
    }else{
      curr_chl <- curr_chl_file[c(181:360,1:180),,,]*1e6*chl_correct ## Extract surface, convert from kg m-3 to mg m-3
    }
    
    surf_chl <- curr_chl[,,1,] # Extract surface values
    curr_chl <- aperm(replicate(dim(curr_chl)[3],surf_chl), c(1,2,4,3)) # Overwrite all non-surface values with surface values (lazy)
    curr_chl <- curr_chl[,,depth_idx,]
    rm(curr_chl_file)
    
    ## NO3
    no3_correct <- esm_scalars[4,i+1]
    curr_no3_file <- var.get.nc(open.nc(list.files(path = paste("./esm_forcings/", curr_esm, "/", sep = ""), pattern = glob2rx(paste("no3*", curr_clim, "*",sep = "")), full.names = TRUE)), "no3")
    curr_no3 <- curr_no3_file[c(181:360,1:180),,depth_idx,]*1e3*no3_correct # to convert from mol m-3 to umol kg-1 (assume 1000kg in 1 m-3)
    rm(curr_no3_file)
    
    ## AOU
    aou_correct <- esm_scalars[6,i+1]
    curr_aou_file <- readRDS(list.files(path = paste("./esm_forcings/", curr_esm, "/", sep = ""), pattern = glob2rx(paste("aou*", curr_clim, "*",sep = "")), full.names = TRUE))
    curr_aou <- curr_aou_file[c(181:360,1:180),,depth_idx,] + aou_correct
    rm(curr_aou_file)
    
    
    curr_data <- expand.grid("Lon" = -179.5:179.5, "Lat" = -89.5:89.5, "logdepth" = log10(curr_depths), "Year" = start_year[j]:end_year[j])
    curr_data <- data.frame("Lon" = curr_data$Lon, "Lat" = curr_data$Lat, "logdepth" = curr_data$logdepth, "depth_m" = 10^curr_data$logdepth, "Year" = curr_data$Year, 
                            "temperature" = as.vector(curr_thetao), "Temp_C" = as.vector(curr_thetao), "SST" = as.vector(SST),
                            "log_nitrate" = log10(as.vector(curr_no3)), "phyc" = as.vector(phyc), "intpp" = as.vector(intpp),
                            "logchlo" = log10(as.vector(curr_chl)), "aou" = as.vector(curr_aou), "BiomassMethod" = "Carbon")
    
    
    curr_data <- curr_data %>% tidyr::drop_na() 
    
    rm(curr_aou, curr_chl, curr_no3, curr_thetao, SST, phyc, intpp)
    
    curr_data <- left_join(curr_data, bathy, by = c("Lon" = "Long", "Lat" = "Lat"))
    
    curr_data <- curr_data %>% dplyr::select(Lon, Lat, Year, depth_m, logdepth, temperature, log_nitrate, logchlo, phyc, intpp,
                                             aou,  SST, Temp_C, BiomassMethod, BATHY)  
    
    years <- unique(curr_data$Year)
    
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = length(years), # Maximum value of the progress bar
                         style = 3)
    
    for(k in 1:length(years)){ ## Vector limit reached when calculations attempted all at once, need to do by year, then bring back together again
      setTxtProgressBar(pb, k)
      year_data1 <- curr_data %>% dplyr::filter(Year == years[k]) %>% group_by(Lon,Lat, Year) %>% 
        dplyr::mutate(depth_m = if_else(depth_m < max(depth_m), depth_m, 200),
                      depth_m = if_else(depth_m > BATHY, BATHY, depth_m),
                      depth_abs = depth_m - dplyr::lag(depth_m, default=0)) 
      
      year_data1$bac_abund_ml <- 10^predict(gm1, year_data1) # Bacteria abundance
      year_data1$bac_cell_carbon_fgC <- predict(glmm1, year_data1, type = "response", re.form=NA) # Bacteria cell carbon
      year_data1$sgr <- 10^(-1.665 + 0.0311771*year_data1$Temp_C + 0.134835*year_data1$logchlo)
      
      microzoo_mgC_m3 <- 10^predict(microzoo_model, year_data1) # Microzoo biomass
      rhizaria_mgC_m3 <- 10^predict(rhizaria_model, year_data1) # Rhizaria biomass
      mesozoo_mgC_m3 <- 10^predict(mesozoo_model, year_data1) # Mesozoo biomass
      macrozoo_mgC_m3 <- 10^predict(macrozoo_model, year_data1) # Macrozoo biomass
      
      year_data1$microzoo_mgC_m3 <- microzoo_mgC_m3
      year_data1$mesozoo_mgC_m3 <- mesozoo_mgC_m3 + rhizaria_mgC_m3
      year_data1$macrozoo_mgC_m3 <- macrozoo_mgC_m3
      year_data1$zooplankton_mgC_m3 <- microzoo_mgC_m3 + rhizaria_mgC_m3 + mesozoo_mgC_m3 + macrozoo_mgC_m3
      rm("microzoo_mgC_m3", "rhizaria_mgC_m3", "mesozoo_mgC_m3", "macrozoo_mgC_m3")
      
      year_data1 <- year_data1 %>% dplyr::mutate(bac_biomass_gCm3 = bac_abund_ml*1e6*bac_cell_carbon_fgC*1e-15,
                                                 bac_growth_gCm3_d = bac_biomass_gCm3*sgr,
                                                 bac_resp_gC_m3_d = bac_growth_gCm3_d*(1-0.14)/0.14,
                                                 bac_resp_urrutia_gC_m3_d = bac_abund_ml*1e6*3.21e11*exp(-0.589/(8.62e-5*(temperature + 273.15)))*1e-15,
                                                 bac_total_gC_m3_d = bac_growth_gCm3_d + bac_resp_gC_m3_d,
                                                 picophy = (0.13*(1-exp(-0.8/0.13*(10^logchlo))))/(10^logchlo),
                                                 nanophy = (0.77*(1-exp(-0.94/0.77*(10^logchlo))))/(10^logchlo) - picophy,
                                                 microphy = ((10^logchlo) - 0.77*(1-exp(-0.94/0.77*(10^logchlo))))/(10^logchlo)) 
      
      year_data1 <- year_data1 %>% group_by(Lon,Lat, Year) %>% 
        dplyr::summarise(bac_biomass_gC_m2 = sum(bac_biomass_gCm3*depth_abs, na.rm = TRUE),
                         bac_cell_carbon_fgC = sum(bac_cell_carbon_fgC*depth_abs, na.rm = TRUE)/sum(depth_abs, na.rm = TRUE),
                         bac_growth_gC_m2_d = sum(bac_growth_gCm3_d*depth_abs, na.rm = TRUE),
                         bac_resp_gC_m2_d = sum(bac_resp_gC_m3_d*depth_abs, na.rm = TRUE),
                         bac_resp_urrutia_gC_m2_d = sum(bac_resp_urrutia_gC_m3_d*depth_abs, na.rm = TRUE),
                         bac_total_gC_m2_d = sum(bac_total_gC_m3_d*depth_abs, na.rm = TRUE),
                         temperature = sum(temperature*depth_abs, na.rm = TRUE)/sum(depth_abs, na.rm = TRUE),
                         nitrate = sum(10^log_nitrate*depth_abs, na.rm = TRUE)/sum(depth_abs, na.rm = TRUE),
                         chlo = sum(10^logchlo*depth_abs, na.rm = TRUE)/sum(depth_abs, na.rm = TRUE),
                         intpp_gC_m2_d = sum(intpp*depth_abs, na.rm = TRUE)/sum(depth_abs, na.rm = TRUE)*12.01*60*60*24, # Convert from moles carbon m2 second to grams carbon m2 day
                         picophyc_vint_gC_m2 = sum(phyc*depth_abs*picophy, na.rm = TRUE)*12.01, # Convert from moles carbon to grams carbon
                         nanophyc_vint_gC_m2 = sum(phyc*depth_abs*nanophy, na.rm = TRUE)*12.01, # Convert from moles carbon to grams carbon
                         microphyc_vint_gC_m2 = sum(phyc*depth_abs*microphy, na.rm = TRUE)*12.01, # Convert from moles carbon to grams carbon
                         phyc_vint_gC_m2 = sum(phyc*depth_abs, na.rm = TRUE)*12.01, # Convert from moles carbon to grams carbon
                         aou = sum(aou*depth_abs, na.rm = TRUE),
                         microzoo_gC_m2 = sum(microzoo_mgC_m3*depth_abs/1e3, na.rm = TRUE),
                         mesozoo_gC_m2 = sum(mesozoo_mgC_m3*depth_abs/1e3, na.rm = TRUE),
                         macrozoo_gC_m2 = sum(macrozoo_mgC_m3*depth_abs/1e3, na.rm = TRUE),
                         zooplankton_gC_m2 = sum(zooplankton_mgC_m3*depth_abs/1e3, na.rm = TRUE), .groups = "drop")
      
      if(k == 1){
        year_data_all1 <- year_data1
      }
      
      if(k != 1){
        year_data_all1 <- rbind(year_data_all1, year_data1)
      }
      
    }
    
    curr_data <- year_data_all1
    rm(year_data_all1)
    
    curr_data <- left_join(curr_data, bathy, by = c("Lon" = "Long", "Lat" = "Lat"))
    curr_data$Experiment <- curr_clim
    curr_data$Model <- curr_esm
    
    if(curr_clim == "historical"){
      spatial_data <- curr_data %>% ungroup() %>% group_by(Lon, Lat, Experiment, Model) %>% dplyr::filter(between(Year,1980,2000)) %>%
        dplyr::summarise(bac_biomass_gC_m2 = mean(bac_biomass_gC_m2, na.rm = TRUE),
                         bac_cell_carbon_fgC = mean(bac_cell_carbon_fgC, na.rm = TRUE),
                         bac_resp_gC_m2_d = mean(bac_resp_gC_m2_d, na.rm = TRUE),
                         bac_resp_urrutia_gC_m2_d = mean(bac_resp_urrutia_gC_m2_d, na.rm = TRUE),
                         bac_growth_gC_m2_d = mean(bac_growth_gC_m2_d, na.rm = TRUE),
                         bac_total_gC_m2_d = mean(bac_total_gC_m2_d, na.rm = TRUE),
                         temperature = mean(temperature, na.rm = TRUE),
                         nitrate = mean(nitrate, na.rm = TRUE),
                         chlo = mean(chlo, na.rm = TRUE),
                         intpp_gC_m2_d = mean(intpp_gC_m2_d, na.rm = TRUE), 
                         phyc_vint_gC_m2 = mean(phyc_vint_gC_m2, na.rm = TRUE),
                         aou = mean(aou, na.rm = TRUE),
                         zooplankton_gC_m2 = mean(zooplankton_gC_m2), .groups = "drop")
    }
    
    if(curr_clim != "historical"){
      spatial_data <- curr_data %>% ungroup() %>% group_by(Lon, Lat, Experiment, Model) %>% dplyr::filter(between(Year,2080,2100)) %>%
        dplyr::summarise(bac_biomass_gC_m2 = mean(bac_biomass_gC_m2, na.rm = TRUE),
                         bac_cell_carbon_fgC = mean(bac_cell_carbon_fgC, na.rm = TRUE),
                         bac_resp_gC_m2_d = mean(bac_resp_gC_m2_d, na.rm = TRUE),
                         bac_resp_urrutia_gC_m2_d = mean(bac_resp_urrutia_gC_m2_d, na.rm = TRUE),
                         bac_growth_gC_m2_d = mean(bac_growth_gC_m2_d, na.rm = TRUE),
                         bac_total_gC_m2_d = mean(bac_total_gC_m2_d, na.rm = TRUE),
                         temperature = mean(temperature, na.rm = TRUE),
                         nitrate = mean(nitrate, na.rm = TRUE),
                         chlo = mean(chlo, na.rm = TRUE),
                         intpp_gC_m2_d = mean(intpp_gC_m2_d, na.rm = TRUE), 
                         phyc_vint_gC_m2 = mean(phyc_vint_gC_m2, na.rm = TRUE),
                         aou = mean(aou, na.rm = TRUE),
                         zooplankton_gC_m2 = mean(zooplankton_gC_m2), .groups = "drop")
    }
    
    
    year_data <- curr_data %>% ungroup() %>% group_by(Year, Experiment, Model) %>% dplyr::filter(Year >= 1980) %>%
      dplyr::summarise(bac_biomass_gC = sum(bac_biomass_gC_m2*area_km2*1e6),
                       bac_cell_carbon_fgC = sum(bac_cell_carbon_fgC*area_km2, na.rm = TRUE)/sum(area_km2*temperature/temperature, na.rm = TRUE),
                       bac_resp_gC_d = sum(bac_resp_gC_m2_d*area_km2*1e6, na.rm = TRUE),
                       bac_resp_urrutia_gC_d = sum(bac_resp_urrutia_gC_m2_d*area_km2*1e6, na.rm = TRUE),
                       bac_growth_gC_d = sum(bac_growth_gC_m2_d*area_km2*1e6, na.rm = TRUE),
                       bac_total_gC_d = sum(bac_total_gC_m2_d*area_km2*1e6, na.rm = TRUE),
                       temperature = sum(temperature*area_km2, na.rm = TRUE)/sum(area_km2*temperature/temperature, na.rm = TRUE),
                       nitrate = sum(nitrate*area_km2, na.rm = TRUE)/sum(area_km2*temperature/temperature, na.rm = TRUE),
                       chlo = sum(chlo*area_km2, na.rm = TRUE)/sum(area_km2*temperature/temperature, na.rm = TRUE),
                       intpp_gC_d = sum(intpp_gC_m2_d*area_km2, na.rm = TRUE)/sum(area_km2*temperature/temperature, na.rm = TRUE), 
                       picophyc_vint_gC = sum(picophyc_vint_gC_m2*area_km2*1e6, na.rm = TRUE),
                       nanophyc_vint_gC = sum(nanophyc_vint_gC_m2*area_km2*1e6, na.rm = TRUE),
                       microphyc_vint_gC = sum(microphyc_vint_gC_m2*area_km2*1e6, na.rm = TRUE),
                       phyc_vint_gC = sum(phyc_vint_gC_m2*area_km2*1e6, na.rm = TRUE),
                       aou = sum(aou*area_km2, na.rm = TRUE)/sum(area_km2*temperature/temperature, na.rm = TRUE),
                       zooplankton_gC = sum(zooplankton_gC_m2*area_km2*1e6, na.rm = TRUE),
                       microzoo_gC = sum(microzoo_gC_m2*area_km2*1e6, na.rm = TRUE),
                       mesozoo_gC = sum(mesozoo_gC_m2*area_km2*1e6, na.rm = TRUE),
                       macrozoo_gC = sum(macrozoo_gC_m2*area_km2*1e6, na.rm = TRUE), .groups = "drop") 
    
    if(curr_clim == "historical"){ # Get mean temperature in 1980-2000 (for fish biomass and respiration deltas)
      start_temp <- year_data %>% ungroup() %>% dplyr::filter(Year <= 2000) %>% dplyr::summarise(temperature = mean(temperature))
    }
    
    year_data <- year_data %>% dplyr::mutate(fish_biomass_lotze = ifelse(temperature >= start_temp$temperature, 5.53*(1-(0.05*(temperature - start_temp$temperature))),
                                                                         5.53*(1+(0.05*(start_temp$temperature - temperature)))),  # 5.53 Gt fish biomass in top 200m from 10.1126/sciadv.abh3732
                                             fish_resp = fish_biomass_lotze*2^(temperature/10),
                                             zoo_resp = zooplankton_gC*2^(temperature/10))
    
    rm(curr_data)
    
    if(j == 1){
      year_data_frame <- year_data
      spatial_data_frame <- spatial_data
    }
    
    if(j != 1){
      year_data_frame <- rbind(year_data_frame, year_data)
      spatial_data_frame <- rbind(spatial_data_frame, spatial_data)
    }
  }
  
  data.table::fwrite(year_data_frame, file = paste("./esm_forcings/processed/", curr_esm, "_year_data_processed_epipelagic.csv", sep = ""), row.names = FALSE)
  data.table::fwrite(spatial_data_frame, file = paste("./esm_forcings/processed/", curr_esm,"_spatial_data_processed_epipelagic.csv", sep = ""), row.names = FALSE)
  
  rm("year_data_frame", "spatial_data_frame")
} 


## Step 2, append all esm data together
esms <- c("GFDL-ESM4", "IPSL-CM6A-LR", "CMCC-ESM2", "MPI-ESM1-2-LR")
year_data_files <- list.files(path = "./esm_forcings/processed", pattern = "year_data", full.names = TRUE)
spatial_data_files <- list.files(path = "./esm_forcings/processed", pattern = "spatial_data", full.names = TRUE)

for(i in 1:length(year_data_files)){
  year_data_frame <- data.table::fread(year_data_files[i])
  spatial_data_frame <- data.table::fread(spatial_data_files[i])
  
  if(i == 1){
    year_data_all <- year_data_frame
    spatial_data_all <- spatial_data_frame
  }
  
  if(i != 1){
    year_data_all <- rbind(year_data_all, year_data_frame)
    spatial_data_all <- rbind(spatial_data_all, spatial_data_frame)
  }
}

data.table::fwrite(year_data_all, file = "./esm_forcings/year_data_processed_epipelagic.csv", row.names = FALSE)
data.table::fwrite(spatial_data_all, file = "./esm_forcings/spatial_data_processed_epipelagic.csv", row.names = FALSE)


