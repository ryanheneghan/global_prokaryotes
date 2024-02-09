## Extended Data Figure 3: Map of prokaryote biomass/all heterotrophs in top 200m

rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)

setwd("~/Desktop/Papers/Bacteria_Census/")

source('./scripts/clean_scripts/plot_map_func.R')

## Import data
bac_epi <- data.table::fread("bacteria_epipelagic.csv")

microzoo_biom <- read.csv("./hatton_data/microzoo_predictions.csv")
macrozoo_biom <- read.csv("./hatton_data/macrozoo_predictions.csv")
mesozoo_biom <- read.csv("./hatton_data/mesozoo_predictions.csv")
rhiz_biom <- read.csv("./hatton_data/rhizaria_predictions.csv")
fish_biom <- read.csv("./hatton_data/fish_predictions.csv")
fish_biom$lat <- -fish_biom$lat

all_heteros <- microzoo_biom %>% dplyr::select(lon, lat, microzoo_biom_gm2_top200m) %>% 
                left_join(., macrozoo_biom %>% dplyr::select(lon, lat, macrozoo_biom_gm2_top200m), by = c("lon" = "lon", "lat" = "lat")) %>%
                left_join(., mesozoo_biom %>% dplyr::select(lon, lat, mesozoo_biom_gm2_top200m), by = c("lon" = "lon", "lat" = "lat")) %>%
                left_join(., rhiz_biom %>% dplyr::select(lon, lat, rhiz_biom_gm2_top200m), by = c("lon" = "lon", "lat" = "lat")) %>%
                left_join(., fish_biom %>% dplyr::select(lon, lat, fish_biom_gm2_1mg_to_1e5g), by = c("lon" = "lon", "lat" = "lat")) %>%
                left_join(., bac_epi %>% dplyr::select(Long, Lat, bac_biomass_gC_m2), by = c("lon" = "Long", "lat" = "Lat")) %>%
                dplyr::rename("microzoo_gm2" = "microzoo_biom_gm2_top200m",
                              "mesozoo_gm2" = "mesozoo_biom_gm2_top200m",
                              "macrozoo_gm2" = "macrozoo_biom_gm2_top200m",
                              "rhiz_gm2" = "rhiz_biom_gm2_top200m",
                              "fish_gm2" = "fish_biom_gm2_1mg_to_1e5g") %>%
              dplyr::mutate(bac_frac = bac_biomass_gC_m2*10/(microzoo_gm2+ mesozoo_gm2+ macrozoo_gm2+ rhiz_gm2+ fish_gm2+bac_biomass_gC_m2*10))


bac_frac_data <- all_heteros %>% dplyr::select(lon,lat,bac_frac)
lonlat <- expand.grid("Long" = -179.5:179.5, "Lat" = -89.5:89.5)
bac_frac_data <- left_join(lonlat, bac_frac_data, by = c("Long" = "lon", "Lat" = "lat"))


plot_bac <- map_figure(bac_frac_data, curr_pal = "BuGn", legend_name = "", plot_lims = c(0.05, 0.45), unlog_leg = FALSE, dp = 2, var_name = "")  + 
  theme(plot.title = element_text(face = "bold", size = 8, vjust = 1, margin = unit(c(0.1,0,0,0), unit = "cm")),
        plot.margin = unit(c(-1,0,-1,0), "lines"))

ggsave(plot = plot_chl, filename = "figures/FigureE3.jpeg", width = 90, height = 60, units = "mm")

