## Figure 2: Bacteria abundance, cell size and total biomass in epi, meso and bathypelagic

library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)


source('./scripts/plot_map_func.R')

## Import data
epipelagic <- data.table::fread("./data/bacteria_epipelagic.csv")
mesopelagic <- data.table::fread("./data/bacteria_mesopelagic.csv")
bathypelagic<- data.table::fread("./data/bacteria_bathypelagic.csv")

#### FIGURE 2
vars <- c("bacabund_m3", "fgC_cell", "bac_biomass_gC_m3")
pals <- c("YlOrRd", "PuBuGn", "viridis")
#legs <- c(expression(paste("m"^{3}," (10"^{11}, ")  ", sep = "")), "fg C  ", expression("mg C m"^-3  ))
legs <- c("", "", "")

titles <- c(expression(paste("Abundance, m"^{-3}," (10"^{11}, ")", sep = "")),
            expression("Cell Carbon, fg C cell"^{-1}), 
            expression(paste("Biomass, mg C m"^-3, sep = "")))
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
                          theme(plot.title = element_text(face = "plain", size = 7, vjust = 1, margin = unit(c(0.1,0,0,0), unit = "cm")),
                                plot.margin = unit(c(-1,0,-2,0), "lines")) + 
                              labs(title = titles[i] )
  plot_list[[i-1 + 4]] <- map_figure(curr_mesop, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[2]], var_name = "") +
                             theme(plot.title = element_blank(), plot.margin = unit(c(-1,0,-2,0), "lines")) #plot_limss[[i]][[1]]
  plot_list[[i-1 + 7]] <- map_figure(curr_bathyp, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[3]], var_name = "") +
                             theme(plot.title = element_blank(), plot.margin = unit(c(-2,0,-2,0), "lines"))
  
}
  
tagss <-c("a", "b", "c", "d", "e", "f", "g", "h", " i")

figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 7),
                            ncol = 3, nrow = 3,   vjust = c(4,4,4,3,3,3,2,2,2), hjust = -2, align = "v")
figure2 <- annotate_figure(figure, 
                           left = text_grob("                 >1000m                                       200-1000m                                           <200m           ", 
                                            size = 7, face = "plain", rot = 90, vjust = 1))

ggsave(plot = figure2, filename = "figures/Figure2.pdf", width = 178, height = 120, units = "mm")


### Save Figure 2, 3 and S5 source data
write.csv(epipelagic, "./source_data/Figure2_3_S5_epipelagic.csv", row.names = FALSE)
write.csv(mesopelagic, "./source_data/Figure2_S5_mesopelagic.csv", row.names = FALSE)
write.csv(bathypelagic, "./source_data/Figure2_S5_bathypelagic.csv", row.names = FALSE)
