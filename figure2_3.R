## Figure 3: Bacteria abundance, cell size and total biomass in epi, meso and bathypelagic

rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)

setwd("~/Desktop/Papers/Bacteria_Census/")

source('./scripts/clean_scripts/plot_map_func.R')

## Import data
area_m2 <- as.vector(as.matrix(t(area(raster()))))*1e6

epipelagic <- data.table::fread("bacteria_epipelagic.csv")
epipelagic$area_m2 <- area_m2
mesopelagic <- data.table::fread("bacteria_mesopelagic.csv")
bathypelagic<- data.table::fread("bacteria_bathypelagic.csv")

#### FIGURE 3
vars <- c("bacabund_m3", "fgC_cell", "bac_biomass_gC_m3")
pals <- c("YlOrRd", "PuBuGn", "viridis")
#legs <- c(expression(paste("m"^{3}," (10"^{11}, ")  ", sep = "")), "fg C  ", expression("mg C m"^-3  ))
legs <- c("", "", "")

titles <- c(expression(bold(paste("Abundance, m"^{3}," (10"^{11}, ")", sep = ""))),
            expression(bold("Cell Carbon, fg C cell"^{-1})), 
            expression(bold(paste("Biomass, mg C m"^-3, sep = ""))))
plot_limss <- list(list(c(2.,5.5), c(0.6,2.), c(0.1,0.9)), list(c(6,12), c(6,12), c(6,12)), list(c(2.,4.8), c(0.6,1.8), c(0.2,0.9)))
#plot_limss <- list(list(c(3,10), c(3,25), c(3,35)), list(c(6,12), c(6,12), c(6,12)), list(c(0.2,0.8), c(0.2,1.8), c(0.1,3.5)))
plot_list <- list()

for(i in 1:length(vars)){
  curr_epip <- epipelagic %>% dplyr::group_by(Long, Lat) %>% dplyr::select(vars[i])
  curr_mesop <- mesopelagic %>% dplyr::group_by(Long, Lat) %>% dplyr::select(vars[i])
  curr_bathyp <- bathypelagic %>% dplyr::group_by(Long, Lat) %>% dplyr::select(vars[i])
  
  if(i == 1){
    curr_epip$bacabund_m3 <- curr_epip$bacabund_m3/1e11
    curr_mesop$bacabund_m3 <- curr_mesop$bacabund_m3/1e11
    curr_bathyp$bacabund_m3 <- curr_bathyp$bacabund_m3/1e11
  }
  
  if(i == 3){
    curr_epip$bac_biomass_gC_m3 <- curr_epip$bac_biomass_gC_m3*1000 # Convert grams to mg (don't change name though)
    curr_mesop$bac_biomass_gC_m3 <- curr_mesop$bac_biomass_gC_m3*1000
    curr_bathyp$bac_biomass_gC_m3 <- curr_bathyp$bac_biomass_gC_m3*1000
  }
  
  plot_list[[i-1 + 1]] <- map_figure(curr_epip, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[1]], var_name = "")  + 
                          theme(plot.title = element_text(face = "bold", size = 8, vjust = 1, margin = unit(c(0.1,0,0,0), unit = "cm")),
                                plot.margin = unit(c(-1,0,-2,0), "lines")) + 
                              labs(title = titles[i] )
  plot_list[[i-1 + 4]] <- map_figure(curr_mesop, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[2]], var_name = "") +
                             theme(plot.title = element_blank(), plot.margin = unit(c(-1,0,-2,0), "lines")) #plot_limss[[i]][[1]]
  plot_list[[i-1 + 7]] <- map_figure(curr_bathyp, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[3]], var_name = "") +
                             theme(plot.title = element_blank(), plot.margin = unit(c(-2,0,-2,0), "lines"))
  
}
  
tagss <-c("a", "b", "c", "d", "e", "f", "g", "h", " i")

figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 8),
                            ncol = 3, nrow = 3,   vjust = c(4,4,4,3,3,3,2,2,2), hjust = -2, align = "v")
figure2 <- annotate_figure(figure, 
                           left = text_grob("                >1000m                                200-1000m                                   <200m           ", 
                                            size = 8, face = "bold", rot = 90, vjust = 1))

ggsave(plot = figure2, filename = "figures/Figure3.jpeg", width = 178, height = 120, units = "mm")

### FIGURE 4
curr_epip <- epipelagic %>% dplyr::group_by(Long, Lat) %>% dplyr::select(bac_c_demand_gC_m3_d)
curr_epip$bac_c_demand_gC_m3_d <- log10(curr_epip$bac_c_demand_gC_m3_d*200*365) # Convert m3 to m2 days to years and log10 (don't change name though)

plot_cd <- map_figure(curr_epip, curr_pal = "PuRd", legend_name = "", plot_lims = c(1, 2.30103), unlog_leg = TRUE, dp = 0, var_name = "")  + 
  theme(plot.title = element_text(face = "bold", size = 8, vjust = 1, margin = unit(c(0.1,0,0,0), unit = "cm")),
        plot.margin = unit(c(-1,0,-1,0), "lines"))

ggsave(plot = plot_cd, filename = "figures/Figure4a.jpeg", width = 90, height = 60, units = "mm")

### Supplementary Figure
vars <- c("temperature_C", "nitrate_umol_kg", "aou_umol_kg")
pals <- c("YlOrRd", "Plasma", "TealGrn")
#legs <- c(expression(paste("m"^{3}," (10"^{11}, ")  ", sep = "")), "fg C  ", expression("mg C m"^-3  ))
legs <- c("", "", "")

titles <- c(expression(bold("Temperature, °C")),
            expression(bold(paste("Nitrate, log"[10], "(",mu, "mol kg"^{-1}, ")",sep = ""))), 
            expression(bold(paste("AOU, ", mu, "mol kg"^{-1}, sep = ""))))
plot_limss <- list(list(c(-2,30), c(-2,18), c(-2,10)), list(log10(c(0.1,50)), log10(c(10,50)), log10(c(15,50))), list(c(-5,300), c(-5,300), c(-5,300)))
#plot_limss <- list(list(c(3,10), c(3,25), c(3,35)), list(c(6,12), c(6,12), c(6,12)), list(c(0.2,0.8), c(0.2,1.8), c(0.1,3.5)))
plot_list <- list()

for(i in 1:length(vars)){
  curr_epip <- epipelagic %>% dplyr::group_by(Long, Lat) %>% dplyr::select(vars[i])
  curr_mesop <- mesopelagic %>% dplyr::group_by(Long, Lat) %>% dplyr::select(vars[i])
  curr_bathyp <- bathypelagic %>% dplyr::group_by(Long, Lat) %>% dplyr::select(vars[i])
  
  if(i == 2){
    curr_epip$nitrate_umol_kg <- log10(curr_epip$nitrate_umol_kg)
    curr_mesop$nitrate_umol_kg <- log10(curr_mesop$nitrate_umol_kg)
    curr_bathyp$nitrate_umol_kg <- log10(curr_bathyp$nitrate_umol_kg)
  }
  
  plot_list[[i-1 + 1]] <- map_figure(curr_epip, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[1]], var_name = "")  + 
    theme(plot.title = element_text(face = "bold", size = 8, vjust = 1, margin = unit(c(0.1,0,0,0), unit = "cm")),
          plot.margin = unit(c(-1,0,-2,0), "lines")) + 
    labs(title = titles[i] )
  plot_list[[i-1 + 4]] <- map_figure(curr_mesop, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[2]], var_name = "") +
    theme(plot.title = element_blank(), plot.margin = unit(c(-1,0,-2,0), "lines")) #plot_limss[[i]][[1]]
  plot_list[[i-1 + 7]] <- map_figure(curr_bathyp, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[3]], var_name = "") +
    theme(plot.title = element_blank(), plot.margin = unit(c(-2,0,-2,0), "lines"))
  
}

tagss <-c("A", "B", "C", "D", "E", "F", "G", "H", " I")

figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 8),
                            ncol = 3, nrow = 3,   vjust = c(4,4,4,3,3,3,2,2,2), hjust = -2, align = "v")
figure2 <- annotate_figure(figure, 
                           left = text_grob("                >1000m                                200-1000m                                   <200m           ", 
                                            size = 8, face = "bold", rot = 90, vjust = 1))

ggsave(plot = figure2, filename = "figures/FigureS3.jpeg", width = 178, height = 120, units = "mm")


