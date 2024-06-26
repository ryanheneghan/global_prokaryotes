### Climate change projections (biomass and respiration)

library(dplyr)
library(ggplot2)
library(raster)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)
library(ggthemes)
library(data.table)

year_data <- read.csv("./data/year_data_processed_epipelagic.csv")

### BIOMASS
tot_bac_biomass1 <- 1.7 # Pg wet biomass, from this study
tot_phy_biomass1 <- 5.4 # Pg wet biomass, from Hatton et al., 2021
tot_zoo_biomass1 <- 7.8 # Pg wet biomass, from Hatton et al., 2021
tot_fish_biomass1 <- 5.5 # Pg wet biomass, from Hatton et al., 2021

#### RESPIRATION
tot_zoo_resp <- 12.4 # from Nowicki and Devries, 2021
tot_fish_resp <- 0.94 # from Bianchi et al. 2020
tot_bac_resp <- 22.6 # 18.5 Pg respiration, from data
tot_bac_growth <- 1.6 # 1.6 Pg production, from data
tot_bac_demand <- tot_bac_resp + tot_bac_growth

### CALCULATE DELTA BIOMASS FOR BACTERIA, ZOOPLANKTON AND FISH, YEAR DATA
hist_data <- year_data %>% dplyr::filter(Year %in% 1980:1999) %>% dplyr::group_by(Model) %>% 
                dplyr::summarise(hist_bac_biomass = mean(bac_biomass_gC),
                                 hist_phy_biomass = mean(phyc_vint_gC),
                                 hist_fish_biomass = mean(fish_biomass_lotze),
                                 hist_zoo_biomass = mean(zooplankton_gC),
                                 hist_bac_resp = mean(bac_resp_gC_d),
                                 hist_zoo_resp = mean(zoo_resp),
                                 hist_fish_resp = mean(fish_resp)) %>% dplyr::ungroup()

year_data2 <- left_join(year_data, hist_data, by = c("Model" = "Model")) %>% 
                    dplyr::select(Year, Experiment, Model, bac_biomass_gC, fish_biomass_lotze,
                                  zooplankton_gC, phyc_vint_gC, hist_bac_biomass,  hist_fish_biomass, hist_zoo_biomass, hist_phy_biomass,
                                  hist_fish_resp, hist_zoo_resp, hist_bac_resp, bac_resp_gC_d, zoo_resp, fish_resp) %>%
                    dplyr::mutate(bac_biomass = (bac_biomass_gC - hist_bac_biomass)/hist_bac_biomass*100,
                                  phy_biomass = (phyc_vint_gC - hist_phy_biomass)/hist_phy_biomass*100,
                                  zoo_biomass = (zooplankton_gC - hist_zoo_biomass)/hist_zoo_biomass*100,
                                  fish_biomass = (fish_biomass_lotze - hist_fish_biomass)/hist_fish_biomass*100,
                                  abs_bac_biomass = -tot_bac_biomass1*(bac_biomass/100),
                                  abs_phy_biomass = -tot_phy_biomass1*(phy_biomass/100),
                                  abs_zoo_biomass = -tot_zoo_biomass1*(zoo_biomass/100),
                                  abs_fish_biomass = -tot_fish_biomass1*(fish_biomass/100),
                                  bac_resp = (bac_resp_gC_d - hist_bac_resp)/hist_bac_resp*100,
                                  zoo_resp = (zoo_resp - hist_zoo_resp)/hist_zoo_resp*100,
                                  fish_resp = (fish_resp - hist_fish_resp)/hist_fish_resp*100,
                                  abs_bac_resp = tot_bac_resp*(bac_resp/100),
                                  abs_fish_resp = tot_fish_resp*(fish_resp/100),
                                  abs_zoo_resp = tot_zoo_resp*(zoo_resp/100))

year_data_gg <- year_data2 %>% dplyr::group_by(Year, Experiment) %>% 
                          summarise(mean_bac_biomass = mean(bac_biomass),
                                    mean_fish_biomass = mean(fish_biomass),
                                    mean_zoo_biomass = mean(zoo_biomass),
                                    mean_phy_biomass = mean(phy_biomass),
                                    mean_abs_bac_biomass = mean(abs_bac_biomass),
                                    mean_abs_phy_biomass = mean(abs_phy_biomass),
                                    mean_abs_zoo_biomass = mean(abs_zoo_biomass),
                                    mean_abs_fish_biomass = mean(abs_fish_biomass),
                                    sd_zoo_biomass = sd(zoo_biomass),
                                    sd_phy_biomass = sd(phy_biomass),
                                    sd_bac_biomass = sd(bac_biomass),
                                    sd_fish_biomass = sd(fish_biomass),
                                    bac_biomass_l = mean_bac_biomass - sd_bac_biomass,
                                    bac_biomass_u = mean_bac_biomass + sd_bac_biomass,
                                    phy_biomass_l = mean_phy_biomass - sd_phy_biomass,
                                    phy_biomass_u = mean_phy_biomass + sd_phy_biomass,
                                    fish_biomass_l = mean_fish_biomass - sd_fish_biomass,
                                    fish_biomass_u = mean_fish_biomass + sd_fish_biomass,
                                    zoo_biomass_l = mean_zoo_biomass - sd_zoo_biomass,
                                    zoo_biomass_u = mean_zoo_biomass + sd_zoo_biomass,
                                    mean_bac_resp = mean(bac_resp),
                                    mean_zoo_resp = mean(zoo_resp),
                                    mean_fish_resp = mean(fish_resp),
                                    sd_bac_resp = sd(bac_resp),
                                    sd_zoo_resp = sd(zoo_resp),
                                    sd_fish_resp = sd(fish_resp),
                                    bac_resp_u = mean_bac_resp + sd_bac_resp,
                                    bac_resp_l = mean_bac_resp - sd_bac_resp,
                                    zoo_resp_u = mean_zoo_resp + sd_zoo_resp,
                                    zoo_resp_l = mean_zoo_resp - sd_zoo_resp,
                                    fish_resp_u = mean_fish_resp + sd_fish_resp,
                                    fish_resp_l = mean_fish_resp - sd_fish_resp,
                                    mean_abs_bac_resp = mean(abs_bac_resp),
                                    mean_abs_fish_resp = mean(abs_fish_resp),
                                    mean_abs_zoo_resp = mean(abs_zoo_resp),
                                    .groups = "drop") %>%
                            dplyr::filter(Experiment %in% c("historical", "ssp370")) 


theme_save <- theme(axis.title = element_text(size = 7, face = "plain"),
                    plot.title = element_blank(),
                    axis.text = element_text(size = 7, colour = "black"),
                    plot.margin = unit(c(0.5,0.8,0.1,0.5), "lines"),
                    legend.title = element_blank(),
                    legend.key.width = unit(0.3, "cm"),
                    legend.key.height = unit(0.15, "cm"),
                    legend.text = element_text(size = 7, face = "plain"))

plot_list <- list()

#colors <- c("Prokaryotes" = "#4274FF", "Phytoplankton" = "#ff9900", "Zooplankton" = "#42B12F", "Fish" = "#EB370D")
#colors2 <- c("#4274FF","#ff9900","#42B12F","#EB370D")

colors <- c("Prokaryotes" = "#882255", "Phytoplankton" = "#DDCC77", "Zooplankton" = "#88CCEE", "Fish" = "#332288")
colors2 <- c("#882255","#DDCC77","#88CCEE","#332288")

plot_list[[1]] <- ggplot(year_data_gg) +
  geom_ribbon(aes(x = Year, ymin=bac_biomass_l,ymax=bac_biomass_u, fill = "Prokaryotes"), alpha=0.2) +
  geom_line(aes(x = Year, y = mean_bac_biomass, color = "Prokaryotes"), linewidth = 1) +
  geom_ribbon(aes(x = Year, ymin=phy_biomass_l,ymax=phy_biomass_u, fill = "Phytoplankton"), alpha=0.2) +
  geom_line(aes(x = Year, y = mean_phy_biomass, color = "Phytoplankton"), linewidth = 1) +
    geom_ribbon(aes(x = Year, ymin=fish_biomass_l,ymax=fish_biomass_u, fill = "Fish"), alpha=0.2) +
    geom_line(aes(x = Year, y = mean_fish_biomass, color = "Fish"), linewidth = 1) +
  geom_ribbon(aes(x = Year, ymin=zoo_biomass_l,ymax=zoo_biomass_u, fill = "Zooplankton"), alpha=0.2) +
  geom_line(aes(x = Year, y = mean_zoo_biomass, color = "Zooplankton"), linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth =0.8)+ theme_bw() + theme_save +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous() + xlab("") + ylab("")+
  coord_cartesian(ylim =c(-13,3), xlim = c(1980,2100))+
  labs(title = "", x = "Year", y = expression(paste(Delta,  " Biomass (%)")), color = "")+  
  scale_color_manual(values = colors,breaks = c("Prokaryotes","Phytoplankton", "Zooplankton", "Fish"),aesthetics = c("colour", "fill")) +
  guides(fill = guide_legend(override.aes = list(size = 1, alpha = 1)), color = "none")


stacked_gg <- year_data_gg %>% dplyr::rename("stack1" = "mean_abs_bac_biomass", "stack2" = "mean_abs_phy_biomass", "stack3" = "mean_abs_zoo_biomass", "stack4" = "mean_abs_fish_biomass")%>% dplyr::select(Year, stack1, stack2, stack3, stack4) %>%
                  gather(variable, value, -Year)


plot_list[[2]] <- ggplot(stacked_gg) + geom_area(aes(x=Year, y = value, fill = variable)) + theme_bw() + theme_save +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + xlab("") + ylab("")+
  coord_cartesian(ylim =c(0,1.75), xlim = c(1980,2100))+
  labs(title = "", x = "Year", y = "Biomass loss (Pg)", color = "")+  
  scale_color_manual(values = colors2, labels = c("Prokaryotes", "Phytoplankton", "Zooplankton","Fish"),aesthetics = c("colour", "fill")) +
  guides(fill = guide_legend(override.aes = list(size = 1, alpha = 1)), color = "none")+
  guides(fill = guide_legend(override.aes = list(alpha = 1)), color = "none")



#### RESPIRATION

plot_list[[3]] <- ggplot(year_data_gg) +
  geom_ribbon(aes(x = Year, ymin=bac_resp_l,ymax=bac_resp_u, fill = "Prokaryotes"), alpha=0.2) +
  geom_line(aes(x = Year, y = mean_bac_resp, color = "Prokaryotes"), linewidth = 1) +
  geom_ribbon(aes(x = Year, ymin=fish_resp_l,ymax=fish_resp_u, fill = "Fish"), alpha=0.2) +
  geom_line(aes(x = Year, y = mean_fish_resp, color = "Fish"), linewidth = 1) +
  geom_ribbon(aes(x = Year, ymin=zoo_resp_l,ymax=zoo_resp_u, fill = "Zooplankton"), alpha=0.2) +
  geom_line(aes(x = Year, y = mean_zoo_resp, color = "Zooplankton"), linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth =0.8)+ theme_bw() + theme_save +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous() + xlab("") + ylab("")+
  coord_cartesian(ylim =c(-3,10), xlim = c(1980,2100))+
  labs(title = "", x = "Year", y = expression(paste(Delta,  " Respiration (%)")), color = "")+  
  scale_color_manual(values = colors,breaks = c("Prokaryotes", "Zooplankton", "Fish"),aesthetics = c("colour", "fill")) +
  guides(fill = guide_legend(override.aes = list(size = 1, alpha = 1)), color = "none")


stacked_gg <- year_data_gg %>% dplyr::rename("stack1" = "mean_abs_bac_resp", "stack2" = "mean_abs_zoo_resp", "stack3" = "mean_abs_fish_resp") %>% dplyr::select(Year, stack1, stack2, stack3) %>%
  gather(variable, value, -Year) 

plot_list[[4]] <- ggplot(stacked_gg) + geom_area(aes(x=Year, y = value, fill = variable)) + theme_bw() + theme_save +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + xlab("") + ylab("")+
  coord_cartesian(ylim =c(0,2.3), xlim = c(1980,2100))+
  labs(title = "", x = "Year", y = expression(paste("Additional respiration (Pg C yr"^{-1}, ")", sep = "")), color = "")+  
  scale_color_manual(values = colors2[c(1,3,4)], labels = c("Prokaryotes", "Zooplankton", "Fish"),aesthetics = c("colour", "fill")) +
  guides(fill = guide_legend(override.aes = list(size = 1, alpha = 1)), color = "none")+
  guides(fill = guide_legend(override.aes = list(alpha = 1)), color = "none")

figure <- ggpubr::ggarrange(plotlist = plot_list, labels = c("a", "b", "c", "d"), font.label = list(size = 7), hjust = -.5, vjust =2, align = 'hv', common.legend = TRUE,legend = "bottom", nrow = 2, ncol = 2)
ggsave(plot = figure, filename = "figures/Figure4.pdf", width = 180, height = 150, units = "mm")


### Save Figure 4 source data
write.csv(year_data_gg, "./source_data/Figure4.csv", row.names = FALSE)

