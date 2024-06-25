### FIGURE S12: Uncertainty distribution of total epipelagic bacteria carbon demand

library(dplyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(raster)
library(ggpubr)

lonlat <- expand.grid("Lon" = -179.5:179.5, "Lat" = -89.5:89.5)
grid_area <- data.frame("Lon" = lonlat$Lon, "Lat" = lonlat$Lat, "area_m2" = 1e6*as.vector(t(as.matrix(area(raster())))))

## Import data
bac_abundance <- read.csv("./data/contemporary_bacteria.csv")
bathy_data <- read.csv("./data/prediction_data.csv")
bathy_data <- bathy_data %>% dplyr::select(Long, Lat, BATHY)

bac_abundance <- left_join(bac_abundance, bathy_data, by = c("Long" = "Long", "Lat" = "Lat"))

euph_data <- bac_abundance %>% filter(depth <= 200) %>% group_by(Long, Lat) %>% mutate(depth = if_else(depth < max(depth), depth, 200), # Stretch deepest interval to go to 200m
                                                               depth = ifelse(depth > BATHY, BATHY, depth), # If deepest interval bottom exceeds bathy, swap to bathy
                                                               depth_abs = depth - dplyr::lag(depth, default=0)) %>%
  dplyr::select(Long, Lat, depth, depth_abs, temperature, chlo_mgm3, bac_biomass_gCm3) %>%
  left_join(grid_area, by = c("Long" = "Lon", "Lat" = "Lat")) %>% dplyr::mutate(sgr_d_mean = 10^(-1.664925+0.031177*temperature+0.134835*log10(chlo_mgm3)))
    
bge_data <- read.csv("./data/robinson_respiration_clean.csv")

n <- 11000
int_coeff <- rnorm(n, mean = -1.664925, sd = 0.118894)

## Remove coefficients outside 95% confidence interval for linear model
int_coeff <- int_coeff[-c(which(int_coeff < -1.664925 - 1.96*0.118894), 
                          which(int_coeff > -1.664925 + 1.96*0.118894))]

n <- 10000
all_bge <- runif(n, min = 0.06, max = 0.27) # Sample from BGE IQR range, assuming uniform distribution
#

store_outputs <- matrix(NA, ncol = n, nrow = 3)  


for(i in 1:n){
  print(i)
 curr_dat <- euph_data %>% dplyr::mutate(sgr_d = 10^(int_coeff[i]+0.031177*temperature+0.134835*log10(chlo_mgm3)), # SPR linear model, from 95% CI for linear model
                                         total_resp_gC_m3_d1 = sgr_d_mean*bac_biomass_gCm3*(1-all_bge[i])/all_bge[i], # Only BGE uncertainty
                                        total_resp_gC_m3_d2 = sgr_d*bac_biomass_gCm3*(1-0.14)/0.14, # Only SPR uncertainty
                                        total_resp_gC_m3_d3 = sgr_d*bac_biomass_gCm3*(1-all_bge[i])/all_bge[i])  %>%  # BGE and SPR uncertainty
                            group_by(Long,Lat, area_m2) %>% dplyr::summarise(bac_resp2_gC_m2_d1 = sum(total_resp_gC_m3_d1*depth_abs, na.rm = TRUE),
                                                                             bac_resp2_gC_m2_d2 = sum(total_resp_gC_m3_d2*depth_abs, na.rm = TRUE),
                                                                             bac_resp2_gC_m2_d3 = sum(total_resp_gC_m3_d3*depth_abs, na.rm = TRUE),
                                                                             .groups = "drop") 

 store_outputs[1,i] <- sum(curr_dat$area_m2*curr_dat$bac_resp2_gC_m2_d1, na.rm = TRUE)*365/1e15 # Only BGE uncertainty
 store_outputs[2,i] <- sum(curr_dat$area_m2*curr_dat$bac_resp2_gC_m2_d2, na.rm = TRUE)*365/1e15 # All SPR uncertainty
 store_outputs[3,i] <- sum(curr_dat$area_m2*curr_dat$bac_resp2_gC_m2_d3, na.rm = TRUE)*365/1e15 # BGE and SPR uncertainty

}

write.csv(store_outputs, "./store_outputs/FigureS12/FigureS12.csv", row.names = FALSE)


### PLOT HISTOGRAMS OF RESULTS
#store_outputs <- as.matrix(read.csv("ed9.csv"))

plot_list <- list()
plot_names <- c("PGE uncertainty", "SPR uncertainty", "PGE + SPR uncertainty")

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", linewidth=0.7),
                         axis.text = element_text(size = 16, face = "bold", colour = "black"),
                         axis.title = element_text(size = 16, face = "bold", colour = "black"),
                         plot.subtitle = element_text(size = 18, face = "bold", colour = "black"),
                         plot.margin = unit(c(0.6,0.3,0.3,0.3), "cm"),
                         legend.key.width = unit(1, "cm"),
                         legend.key.height = unit(2, "cm"),
                         legend.text = element_text(size = 20, face = "bold"),
                         legend.title = element_text(size = 20, face = "bold"),
                         legend.position = 'bottom'))

for(i in 1:length(plot_names)){
  gg_data <- data.frame("dd" = store_outputs[i,])
  
  plot_list[[i]] <- ggplot(data = gg_data, aes(x=dd)) + geom_histogram(aes(y=after_stat(count / sum(count))), bins = 50, col = "black", fill = "lightblue") + theme_bw() +
    theme_opts+ scale_x_continuous(expand = c(0,0), limits = c(0,105)) + scale_y_continuous(expand = c(0,0), limits = c(0,0.155)) +
    labs(x = expression(bold(paste("Prokaryotic Respiration (Pg C yr"^{-1}, ")", sep = ""))),
         y = "Proportion", subtitle = plot_names[i]) +
    geom_vline(aes(xintercept=22.6),
               color="red", linetype="dashed", linewidth=1.2)
  
}

tagss <-c("a", "b", "c")
figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1 , ncol = 3, nrow = 1)
ggsave(plot = figure, filename = "figures/FigureS12.jpeg", width = 470, height = 130, units = "mm")
