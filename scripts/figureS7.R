## SUPPLEMENTARY FIGURE S7, specific-production rate statistical model

library(dplyr)
library(tidyr)
library(viridis)
library(ggpointdensity)
library(ggplot2)
library(ggpubr)
library(lme4)
library(visreg)
library(xlsx)

#### FIGURE OF FINAL PARAMETRIC MODEL FOR SPR
all_data <- read.csv("./data/raw_observations/prokaryote_sgr.csv")
all_data <- all_data %>% dplyr::mutate(across(-c(source), as.numeric))
all_data$logSGR <- log10(all_data$SGR)
all_data$log10CHL <- log10(all_data$CHL)

lmer_temp <- lmer(logSGR ~ SST + (1|source), data = all_data) 
lmer_chl <- lmer(logSGR ~ SST + log10CHL + (1|source), data = all_data) 

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", linewidth=0.7),
                         axis.text = element_text(size = 18, face = "bold", colour = "black"),
                         axis.title = element_text(size = 20, face = "bold", colour = "black"),
                         plot.title = element_text(size = 18, face = "bold.italic", colour = "black"),
                         plot.margin = unit(c(0.6,0.3,0.3,0.3), "cm"),
                         legend.key.width = unit(1, "cm"),
                         legend.key.height = unit(2, "cm"),
                         legend.text = element_text(size = 20, face = "bold"),
                         legend.title = element_text(size = 20, face = "bold"),
                         legend.position = 'bottom'))

plot_list <- list()

tt = visreg(lmer_temp, xvar = "SST", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("Temp_C" = tt$fit$SST, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
figS7a_line <- line_dat %>% rename(log10_SPR_d_mean = "Value", log10_SPR_d_95lowerbound = "Lower", log10_SPR_d_95upperbound = "Upper") 

dot_dat <- data.frame("Temp_C" = tt$res$SST, "Value" = tt$res$visregRes)
figS7a_points <- dot_dat %>% rename(log10_SPR_d = "Value")

x <- densCols(dot_dat$Temp_C,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[1]] <- ggplot() + geom_point(data=dot_dat, aes(x=Temp_C, y = Value, colour = dens), size = 1)+ 
  geom_line(data=line_dat, aes(x=Temp_C, y = Value), col = "red", linewidth = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = Temp_C, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,80,160,250), labels = c("1", "80",  "160",  "≥250"))+
  xlab(expression(bold("Temperature, °C" )))+
  ylab(expression(bold(paste("log"[10], "(SPR, d"^{-1}, ")", sep = "")))) + 
  theme_opts

tt = visreg(lmer_chl, xvar = "log10CHL", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("log10CHL" = tt$fit$log10CHL, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
figS7b_line <- line_dat %>% rename(log10_SPR_d_mean = "Value", log10_SPR_d_95lowerbound = "Lower", log10_SPR_d_95upperbound = "Upper") 

dot_dat <- data.frame("log10CHL" = tt$res$log10CHL, "Value" = tt$res$visregRes)
figS7b_points <- dot_dat %>% rename(log10_SPR_d = "Value")

x <- densCols(dot_dat$log10CHL,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[2]] <- ggplot() + geom_point(data=dot_dat, aes(x=log10CHL, y = Value, colour = dens), size = 1)+ 
  geom_line(data=line_dat, aes(x=log10CHL, y = Value), col = "red", linewidth = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = log10CHL, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,80,160,250), labels = c("1", "80",  "160",  "≥250"))+
  xlab(expression(bold(paste("log"[10], "(Chlorophyll, mg m"^{-3}, ")", sep = ""))))+
  ylab(expression(bold(paste("log"[10], "(SPR, d"^{-1}, ")", sep = "")))) + 
  theme_opts


tagss <-c("a", "b")
figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1 , ncol = 2, nrow = 1, legend = "right", common.legend = TRUE)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS7.png", sep = ""), width = 1000, height = 400)


### Save Figure S7 source data
wb <- createWorkbook(type = "xls")

sh2 <- createSheet(wb,sheetName = "Figure S7a_Points")
addDataFrame(data.frame("Figure S7a_Points"=c("Figure S7a_Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(figS7a_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S7a_Line")
addDataFrame(data.frame("Figure S7a_Line"=c("Figure S7a_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(figS7a_line,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S7b_Points")
addDataFrame(data.frame("Figure S7b_Points"=c("Figure S7b_Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(figS7b_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S7b_Line")
addDataFrame(data.frame("Figure S7b_Line"=c("Figure S7b_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(figS7b_line,sheet = sh2,row.names = FALSE,startRow = 2)

saveWorkbook(wb,"./source_data/FigureS7.xls")

