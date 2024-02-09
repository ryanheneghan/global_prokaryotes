### Climate change projections (biomass and respiration)

rm(list=ls())
library(dplyr)
library(ggplot2)
library(raster)
library(dplyr)
library(reshape2)
library(ggpubr)
library(ggthemes)
library(data.table)

setwd("~/Desktop/Papers/Bacteria_Census/")

load("lm_model_final.RData")
load("glmm_model1_carbon.RData")

pred_data <- data.frame("temperature" = 0:30, "Temp_C" = 0:30, 
                        "logdepth" = 0, "log_nitrate" = 0,
                        "aou" = 0, "logchlo" = 0)

pred_data <- pred_data %>% dplyr::mutate(abund = 10^predict(gm1, pred_data),
                                         size = predict(glmm1, pred_data, type = "response", re.form=NA),
                                         biomass = abund*size)

pred_data$biom_norm <- (pred_data$biomass/pred_data$biomass[1] - 1)*100

theme_save <- theme(axis.title = element_text(size = 7, face = "bold"),
                    plot.title = element_blank(),
                    axis.text = element_text(size = 7, colour = "black"),
                    plot.margin = unit(c(1,1.4,0.1,0.5), "lines"),
                    legend.title = element_blank(),
                    legend.key.width = unit(0.3, "cm"),
                    legend.key.height = unit(0.3, "cm"),
                    legend.text = element_text(size = 7, face = "bold"))

plot1 <- ggplot(pred_data) + geom_line(aes(x=temperature, y = biom_norm), linewidth = 0.7) + 
  theme_bw() + theme_save +   geom_hline(yintercept = 0, color = "black", linetype = "dashed", linewidth =0.8)+ theme_bw() + theme_save +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous() + xlab("") + ylab("")+
  coord_cartesian(ylim =c(-13,20), xlim = c(0,30))+
  labs(title = "", x = expression(bold("Temperature, °C" )), y = expression(bold(paste(Delta,  " Biomass (%)"))), color = "")
  
ggsave(plot = plot1, filename = "figures/ExtendedDataFigure2.jpeg", width = 90, height = 75, units = "mm")

