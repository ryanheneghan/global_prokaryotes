## Extended Data Figure 2: Environmental variables in epi, meso and bathypelagic

library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)

source('./scripts/plot_map_func.R')

## Import data
epipelagic <- data.table::fread("./data/bacteria_epipelagic.csv")
mesopelagic <- data.table::fread("./data/bacteria_mesopelagic.csv")
bathypelagic<- data.table::fread("./data/bacteria_bathypelagic.csv")

vars <- c("temperature_C", "nitrate_umol_kg", "aou_umol_kg")
pals <- c("OrRd", "Plasma", "TealGrn")
#legs <- c(expression(paste("m"^{3}," (10"^{11}, ")  ", sep = "")), "fg C  ", expression("mg C m"^-3  ))
legs <- c("", "", "")

titles <- c(expression(bold("Temperature, °C")),
            expression(bold(paste("Nitrate, " ,mu, "mol kg"^{-1}, sep = ""))), 
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
    theme(plot.title = element_blank(), plot.margin = unit(c(-1,0,-2,0), "lines")) 
  plot_list[[i-1 + 7]] <- map_figure(curr_bathyp, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[3]], var_name = "") +
    theme(plot.title = element_blank(), plot.margin = unit(c(-2,0,-2,0), "lines"))
  
  if(i == 2){
  plot_list[[i-1 + 1]] <- map_figure(curr_epip, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[1]], unlog_leg = TRUE, var_name = "")  + 
    theme(plot.title = element_text(face = "bold", size = 8, vjust = 1, margin = unit(c(0.1,0,0,0), unit = "cm")),
          plot.margin = unit(c(-1,0,-2,0), "lines")) + 
    labs(title = titles[i] )
  plot_list[[i-1 + 4]] <- map_figure(curr_mesop, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[2]], unlog_leg = TRUE, var_name = "") +
    theme(plot.title = element_blank(), plot.margin = unit(c(-1,0,-2,0), "lines")) 
  plot_list[[i-1 + 7]] <- map_figure(curr_bathyp, curr_pal = pals[i], legend_name = legs[i], plot_lims = plot_limss[[i]][[3]], unlog_leg = TRUE, var_name = "") +
    theme(plot.title = element_blank(), plot.margin = unit(c(-2,0,-2,0), "lines"))
  }
  
}

tagss <-c("a", "b", "c", "d", "e", "f", "g", "h", " i")

figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 8),
                            ncol = 3, nrow = 3,   vjust = c(4,4,4,3,3,3,2,2,2), hjust = -2, align = "v")
figure2 <- annotate_figure(figure, 
                           left = text_grob("                >1000m                                200-1000m                                   <200m           ", 
                                            size = 8, face = "bold", rot = 90, vjust = 1))

ggsave(plot = figure2, filename = "figures/FigureED2.jpeg", width = 178, height = 120, units = "mm")

