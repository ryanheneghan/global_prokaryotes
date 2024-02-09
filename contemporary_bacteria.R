### Calculate contemporary bacteria abundance, cell size and biomass using WOA and modis measurements

rm(list = ls())

library(dplyr)
library(tidyr)
library(stats)
library(lme4)
library(reshape2)

setwd("~/Desktop/Papers/Bacteria_Census")

load("lm_model_final.RData") # bacteria abundance model
load("glmm_model1_carbon.RData") # bacteria cell carbon model

## LOAD WOA DATA, AND SATELLITE CHLO, BATHY OF EACH CELL
nitrate <- read.csv("./stat_model/nitrate/woa18_all_n00an01.csv", skip = 1, header = TRUE, check.names = FALSE)
aou <- read.csv("./stat_model/aou/woa18_all_A00an01.csv",skip = 1, header = TRUE, check.names = FALSE)
temperature <- read.csv("./stat_model/temperature/woa18_decav81B0_t00an01.csv", skip = 1, header = TRUE, check.names = FALSE)
colnames(nitrate) <- colnames(aou) <- colnames(temperature) <- c("Lat","Long","0",colnames(nitrate)[-c(1:3)])

nitrate <- nitrate %>% reshape2::melt(id.vars = c("Long", "Lat")) %>% rename(depth = variable, nitrate = value) %>% mutate_if(is.factor, ~ as.numeric(as.character(.x)))
aou <- aou %>% reshape2::melt(id.vars = c("Long", "Lat")) %>% rename(depth = variable, aou = value) %>% mutate_if(is.factor, ~ as.numeric(as.character(.x)))
temperature <- temperature %>% reshape2::melt(id.vars = c("Long", "Lat")) %>% rename(depth = variable, temperature = value) %>% mutate_if(is.factor, ~ as.numeric(as.character(.x)))

lonlat = expand.grid("Lon" = -179.5:179.5, "Lat" = -89.5:89.5)
chlo_nc <- var.get.nc(open.nc("./stat_model/chlorophyll/chlo_modis_aqua_average_2002-2019.nc"), "chlor_a")
chlo_nc <- chlo_nc[c(181:360, 1:180),]
chlo_nc <- data.frame("Long" = lonlat$Lon, "Lat" = lonlat$Lat, "chlo_mgm3" = as.vector(chlo_nc))

all_data <- nitrate
all_data$aou <- aou$aou
all_data$temperature <- temperature$temperature

all_data <- left_join(all_data, chlo_nc, by = c("Long" = "Long", "Lat" = "Lat"))
all_data <- all_data %>% drop_na() %>% dplyr::mutate(depth = depth + 0.5,
                                                     log_nitrate = log10(nitrate), 
                                                     logdepth = log10(depth),
                                                     logchlo = log10(chlo_mgm3),
                                                     Temp_C = temperature)

all_data$log_bacabund_ml <- predict(gm1, all_data)
all_data$bac_cell_carbon_fgC <- predict(glmm1, all_data, type = "response", re.form=NA)
all_data <- all_data %>% dplyr::mutate(bac_biomass_gCm3 = 10^log_bacabund_ml*1e6*bac_cell_carbon_fgC*1e-15)

write.csv(all_data, "contemporary_bacteria.csv", row.names = FALSE)

