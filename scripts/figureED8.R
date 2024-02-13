### Extended data figure 8: Climate change projections on temperature and bacteria cell-carbon

library(dplyr)
library(ggplot2)
library(raster)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr)
library(ggthemes)

year_data <- read.csv("./data/year_data_processed_epipelagic.csv")

### CALCULATE DELTA BIOMASS AND RESPIRATION FOR BACTERIA, YEAR DATA
hist_data <- year_data %>% dplyr::filter(Year %in% 1980:1999) %>% dplyr::group_by(Model) %>% 
                dplyr::summarise(hist_bac_size = mean(bac_cell_carbon_fgC),
                                 hist_temp = mean(temperature),
                                 hist_nitrate = mean(nitrate),
                                 hist_chlo = mean(chlo)) %>% dplyr::ungroup()

year_data2 <- left_join(year_data, hist_data, by = c("Model" = "Model")) %>% 
                    dplyr::select(Year, Experiment, Model, bac_cell_carbon_fgC, temperature, 
                    nitrate, chlo, hist_bac_size, hist_temp, hist_nitrate, hist_chlo) %>%
                    dplyr::mutate(bac_size = (bac_cell_carbon_fgC - hist_bac_size)/hist_bac_size*100,
                                  temperature = (temperature - hist_temp),
                                  chlo = (chlo - hist_chlo)/hist_chlo*100,
                                  nitrate = (nitrate - hist_nitrate)/hist_nitrate*100)

year_data_gg <- year_data2 %>% dplyr::group_by(Year, Experiment) %>% 
                          summarise(mean_bac_size = mean(bac_size),
                                    mean_temperature = mean(temperature),
                                    mean_chlo = mean(chlo),
                                    mean_nitrate = mean(nitrate),
                                    sd_bac_size = sd(bac_size),
                                    sd_temperature = sd(temperature),
                                    sd_chlo = sd(chlo),
                                    sd_nitrate = sd(nitrate),
                                    bac_size_l = mean_bac_size - sd_bac_size,
                                    bac_size_u = mean_bac_size + sd_bac_size,
                                    temperature_l = mean_temperature - sd_temperature,
                                    temperature_u = mean_temperature + sd_temperature,
                                    chlo_l = mean_chlo - sd_chlo,
                                    chlo_u = mean_chlo + sd_chlo,
                                    nitrate_l = mean_nitrate - sd_nitrate,
                                    nitrate_u = mean_nitrate + sd_nitrate,
                                    .groups = "drop") %>%
                            dplyr::filter(Experiment %in% c("historical", "ssp370")) 


theme_save <- theme(axis.title = element_text(size = 7, face = "bold"),
                    plot.title = element_blank(),
                    axis.text = element_text(size = 7, colour = "black"),
                    plot.margin = unit(c(1,1.4,0.1,0.5), "lines"),
                    legend.title = element_blank(),
                    legend.key.width = unit(0.3, "cm"),
                    legend.key.height = unit(0.3, "cm"),
                    legend.text = element_text(size = 7, face = "bold"))

plot_list <- list()

plot_list[[1]] <- ggplot(year_data_gg) +
  geom_ribbon(aes(x = Year, ymin=temperature_l,ymax=temperature_u),  fill = "gray", alpha=0.7) +
  geom_line(aes(x = Year, y = mean_temperature), color = "black", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth =0.8)+ theme_bw() + theme_save +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous() + xlab("") + ylab("")+
  coord_cartesian(ylim =c(0,2.2), xlim = c(1980,2100))+
  labs(title = "", x = "Year", y = expression(bold(paste(Delta,  " Temperature, °C"))), color = "")

plot_list[[2]] <- ggplot(year_data_gg) +
  geom_ribbon(aes(x = Year, ymin=bac_size_l,ymax=bac_size_u),  fill = "gray", alpha=0.7) +
  geom_line(aes(x = Year, y = mean_bac_size), color = "black", linewidth = 1)+
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth =0.8)+ theme_bw() + theme_save +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous() + xlab("") + ylab("")+
  coord_cartesian(ylim =c(-4.5,0.2), xlim = c(1980,2100))+
  labs(title = "", x = "Year", y = expression(bold(paste(Delta,  " Cell Carbon (%)"))), color = "")

figure <- ggpubr::ggarrange(plotlist = plot_list, labels = c("a", "b"), font.label = list(size = 8), hjust = -.5, vjust =2, align = 'hv', nrow = 1, ncol = 2)
ggsave(plot = figure, filename = "figures/FigureED8.jpeg", width = 178, height = 75, units = "mm")


