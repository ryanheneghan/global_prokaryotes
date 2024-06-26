## Figure S20: GAM for SPR

library(dplyr)
library(tidyr)
library(viridis)
library(ggpointdensity)
library(ggplot2)
library(ggpubr)
library(visreg)
library(mgcv)
library(xlsx)

df <- read.csv('./data/raw_observations/prokaryote_sgr.csv', header=TRUE)
df <- df %>% dplyr::mutate(log10SGR = log10(SGR),
                           log10CHL = log10(CHL),
                           log10Depth = log10(Depth))


gam_temp <- gam(log10SGR ~ s(SST) + source, data = df)
gam_chl <- gam(log10SGR ~ s(SST) + s(log10CHL) + source, data = df)

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


# PLOT MAIN EFFECTS
plot_list <- list()

tt = visreg(gam_temp, xvar = "SST", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("Temp_C" = tt$fit$SST, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
fig1a_line <- line_dat %>% rename(log10_SPR_d_mean = "Value", log10_SPR_d_95lowerbound = "Lower", log10_SPR_d_95upperbound = "Upper") 

dot_dat <- data.frame("Temp_C" = tt$res$SST, "Value" = tt$res$visregRes)
fig1a_points <- dot_dat %>% rename(log10_SPR_d = "Value")

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

tt = visreg(gam_chl, xvar = "log10CHL", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("log10CHL" = tt$fit$log10CHL, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
fig1b_line <- line_dat %>% rename(log10_SPR_d_mean = "Value", log10_SPR_d_95lowerbound = "Lower", log10_SPR_d_95upperbound = "Upper") 

dot_dat <- data.frame("log10CHL" = tt$res$log10CHL, "Value" = tt$res$visregRes)
fig1b_points <- dot_dat %>% rename(log10_SPR_d = "Value")

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

tt = visreg(gam_temp, xvar = "source", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("Source" = tt$fit$source, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
fig1c_line <- line_dat %>% rename(log10_SPR_d_mean = "Value", log10_SPR_d_95lowerbound = "Lower", log10_SPR_d_95upperbound = "Upper") 

dot_dat <- data.frame("Source" = tt$res$source, "Value" = tt$res$visregRes)
fig1c_points <- dot_dat %>% rename(log10_SPR_d = "Value")

x <- densCols(dot_dat$Source ,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[3]] <- ggplot() + geom_point(data=dot_dat, aes(x=Source, y = Value, colour = dens), size = 1)+ 
  geom_point(data=line_dat, aes(x=Source, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_pointrange(data = line_dat, aes(x=Source, y = Value, ymin=Lower, ymax=Upper), col = "red")+
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,80,160,250), labels = c("1", "80",  "160",  "≥250"))+
  scale_x_discrete("Dataset", breaks = c("Blanes_Bay", "HotMix", "Latitud", "Malaspina", "Kirchman_2009"), labels = c("Blanes Bay", "HotMix", "Latitud", "Malaspina-2010", "Polar")) +
  ylab(expression(bold(paste("log"[10], "(SPR, d"^{-1}, ")", sep = ""))))  +
  theme_opts + theme(axis.text.x = element_text(size = 10, face = "bold", colour = "black", angle = 20, vjust = 0.7))

tagss <-c("a", "b", "c")
figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1 , ncol = 2, nrow = 2, legend = "right", common.legend = TRUE)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS20.png", sep = ""), width = 1000, height = 800)


### Save Figure S20 source data
wb <- createWorkbook(type = "xls")
sh1 <- createSheet(wb,sheetName = "Figure S20a_Line")
addDataFrame(data.frame("Figure S20a_Line"=c("Figure S20a_Line")),sheet = sh1,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1a_line,sheet = sh1,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S20a_Points")
addDataFrame(data.frame("Figure S20a_Points"=c("Figure S20a_Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1a_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S20b_Line")
addDataFrame(data.frame("Figure S20b_Line"=c("Figure S20b_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1b_line,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S20b_Points")
addDataFrame(data.frame("Figure S20b Points"=c("Figure S20b Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1b_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S20c_Line")
addDataFrame(data.frame("Figure S20c_Line"=c("Figure S20c_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1c_line,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S20c_Points")
addDataFrame(data.frame("Figure S20c Points"=c("Figure S20c Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1c_points,sheet = sh2,row.names = FALSE,startRow = 2)

saveWorkbook(wb,"./source_data/FigureS20.xls")
