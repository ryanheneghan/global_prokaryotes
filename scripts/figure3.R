### FIGURE 3: Bacteria carbon demand in epipelagic

library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)

source('./scripts/plot_map_func.R')

## Import data
epipelagic <- data.table::fread("./data/bacteria_epipelagic.csv")

curr_epip <- epipelagic %>% dplyr::group_by(Long, Lat) %>% dplyr::select(bac_c_demand_gC_m3_d)
curr_epip$bac_c_demand_gC_m3_d <- log10(curr_epip$bac_c_demand_gC_m3_d*200*365) # Convert m3 to m2 days to years and log10 (don't change name though)

plot_cd <- map_figure(curr_epip, curr_pal = "PuRd", legend_name = "", plot_lims = c(1, 2.30103), unlog_leg = TRUE, dp = 0, var_name = "")  + 
  theme(plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), "lines"))

ggsave(plot = plot_cd, filename = "figures/Figure3.pdf", width = 88, height = 55, units = "mm")

