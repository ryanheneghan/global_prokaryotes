## Figure S8: Map of surface chlorophyll

library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)

source('./scripts/plot_map_func.R')

## Import data
chl_data <- read.csv("./data/surface_chl.csv")

plot_chl <- map_figure(chl_data, curr_pal = "Emrld", legend_name = "", plot_lims = c(-1.5, 0.69897), unlog_leg = TRUE, dp = 2, var_name = "")  + 
  theme(plot.title = element_text(face = "bold", size = 8, vjust = 1, margin = unit(c(0.1,0,0,0), unit = "cm")),
        plot.margin = unit(c(-1,0,-1,0), "lines"))

ggsave(plot = plot_chl, filename = "figures/FigureS8.jpeg", width = 90, height = 60, units = "mm")

write.csv(chl_data, "./source_data/FigureS8.csv", row.names = FALSE)
