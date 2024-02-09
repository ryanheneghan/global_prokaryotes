rm(list = ls())

library(dplyr)
library(tidyr)
library(data.table)
library(raster)
library(ggplot)

setwd("~/Desktop/Papers/Bacteria_Census")

### PROCESS ESM DATA
esm_data <- read.csv("./esm_forcings/year_data_processed_epipelagic.csv")
esm_data <- esm_data %>% dplyr::select(Year, Experiment, Model, temperature, bac_biomass_gC,
                                       picophyc_vint_gC, nanophyc_vint_gC, microphyc_vint_gC,
                                       microzoo_gC, mesozoo_gC, macrozoo_gC, fish_biomass_lotze)

hist_data <- esm_data %>% dplyr::filter(Year %in% 1980:1999) %>% dplyr::group_by(Model) %>% 
                  summarise_if(is.numeric, mean, na.rm = TRUE) %>% dplyr::select(-Year) %>% dplyr::ungroup()

all_data <- left_join(esm_data, hist_data, by = c("Model" = "Model")) %>%
  dplyr::mutate(temperature = temperature.x - temperature.y,
                bac = (bac_biomass_gC.x - bac_biomass_gC.y)/bac_biomass_gC.y*100,
                picophy = (picophyc_vint_gC.x - picophyc_vint_gC.y)/picophyc_vint_gC.y*100,
                nanophy = (nanophyc_vint_gC.x - nanophyc_vint_gC.y)/nanophyc_vint_gC.y*100,
                microphy = (microphyc_vint_gC.x - microphyc_vint_gC.y)/microphyc_vint_gC.y*100,
                microzoo = (microzoo_gC.x - microzoo_gC.y)/microzoo_gC.y*100,
                mesozoo = (mesozoo_gC.x - mesozoo_gC.y)/mesozoo_gC.y*100,
                macrozoo = (macrozoo_gC.x - macrozoo_gC.y)/macrozoo_gC.y*100,
                fish = (fish_biomass_lotze.x - fish_biomass_lotze.y)/fish_biomass_lotze.y*100) %>%
                dplyr::select(-contains(c('.x', '.y'))) %>% group_by(Year) %>%
                summarise_if(is.numeric, mean, na.rm = TRUE)
            
bac_temp <- 1+lm(bac ~ temperature, data = all_data)$coefficients[2]/100*2
picophy_temp <- 1+lm(picophy ~ temperature, data = all_data)$coefficients[2]/100*2
nanophy_temp <- 1+lm(nanophy ~ temperature, data = all_data)$coefficients[2]/100*2
microphy_temp <- 1+lm(microphy ~ temperature, data = all_data)$coefficients[2]/100*2
microzoo_temp <- 1+lm(microzoo ~ temperature, data = all_data)$coefficients[2]/100*2
mesozoo_temp <- 1+lm(mesozoo ~ temperature, data = all_data)$coefficients[2]/100*2
macrozoo_temp <- 1+lm(macrozoo ~ temperature, data = all_data)$coefficients[2]/100*2
fish_temp <- 1+lm(fish ~ temperature, data = all_data)$coefficients[2]/100*2

hatton_data <- read.csv("./hatton_data/summary_biomass_top200_table_long.csv")
base_ss <- hatton_data %>% dplyr::group_by(log10_Size_Midpoint_g) %>% dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)
temp_ss <- hatton_data %>% dplyr::mutate(Biomass_Pg_wet_weight_estimate = ifelse(Group == "Hetero_Bacteria", Biomass_Pg_wet_weight_estimate*bac_temp, Biomass_Pg_wet_weight_estimate),
                                         Biomass_Pg_wet_weight_estimate = ifelse(Group == "Picophytoplankton", Biomass_Pg_wet_weight_estimate*picophy_temp, Biomass_Pg_wet_weight_estimate),
                                         Biomass_Pg_wet_weight_estimate = ifelse(Group == "Nanophytoplankton", Biomass_Pg_wet_weight_estimate*nanophy_temp, Biomass_Pg_wet_weight_estimate),
                                         Biomass_Pg_wet_weight_estimate = ifelse(Group == "Microphytoplankton", Biomass_Pg_wet_weight_estimate*microphy_temp, Biomass_Pg_wet_weight_estimate),
                                         Biomass_Pg_wet_weight_estimate = ifelse(Group == "Microzooplankton", Biomass_Pg_wet_weight_estimate*microzoo_temp, Biomass_Pg_wet_weight_estimate),
                                         Biomass_Pg_wet_weight_estimate = ifelse(Group == "Mesozooplankton", Biomass_Pg_wet_weight_estimate*mesozoo_temp, Biomass_Pg_wet_weight_estimate),
                                         Biomass_Pg_wet_weight_estimate = ifelse(Group == "Macrozooplankton", Biomass_Pg_wet_weight_estimate*macrozoo_temp, Biomass_Pg_wet_weight_estimate),
                                         Biomass_Pg_wet_weight_estimate = ifelse(Group == "Fish", Biomass_Pg_wet_weight_estimate*fish_temp, Biomass_Pg_wet_weight_estimate),
                                         Biomass_Pg_wet_weight_estimate = ifelse(Group == "Mammals", Biomass_Pg_wet_weight_estimate*fish_temp, Biomass_Pg_wet_weight_estimate)) %>%
                            group_by(log10_Size_Midpoint_g) %>% dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)

lm(log10(base_ss$Biomass_Pg_wet_weight_estimate) ~ base_ss$log10_Size_Midpoint_g)
lm(log10(temp_ss$Biomass_Pg_wet_weight_estimate) ~ temp_ss$log10_Size_Midpoint_g)

plot(base_ss$log10_Size_Midpoint_g, log10(base_ss$Biomass_Pg_wet_weight_estimate), type = "l", xlab = "log10(Body Size, g)", ylab = "log10(Biomass, Pg)", main = "Blue = 2C warming, \n Black = No warming")
lines(temp_ss$log10_Size_Midpoint_g, log10(temp_ss$Biomass_Pg_wet_weight_estimate), col = "blue")

base_ss$set <- "Original"
temp_ss$set <- "2C"

all_ss <- rbind(base_ss, temp_ss)

ggplot(data = all_ss, aes(x = log10_Size_Midpoint_g, y = Biomass_Pg_wet_weight_estimate)) + 
  geom_bar(stat = "identity", fill = "grey", colour = "red")+
  geom_bar(data = all_ss, aes(x = log10_Size_Midpoint_g, y = set),
           stat = "identity", fill = "grey", colour = "blue")
