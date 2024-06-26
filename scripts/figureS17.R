### Figure S17: GAM MODEL FOR PROK CELL CARBON 

library(dplyr)
library(tidyr)
library(viridis)
library(ggpointdensity)
library(ggplot2)
library(ggpubr)
library(mgcv)
library(visreg)
library(xlsx)

df <- read.csv('./data/raw_observations/prokaryote_cellcarbon.csv', header=TRUE)
df <- df %>% mutate(Depth = as.numeric(Depth),
                    log10_depth = log10(Depth),
                    log10_BT = log10(as.numeric(All_BT)),
                    log10_fg = log10(fgC_cell),
                    log10_BT2 = log10_BT^2,
                    log10_BT3 = log10_BT^3,
                    Station = ifelse(nchar(as.character(Station)) == 1, paste("00", as.character(Station), sep = ""), 
                                     ifelse(nchar(as.character(Station)) == 2, paste("0", as.character(Station), sep = ""), as.character(Station))),
                    Leg = as.character(Leg),
                    Leg_5 = ifelse(Leg == "5", "Yes", "No")) %>% 
  drop_na(Leg, Temp_C, log10_depth, log10_BT)#%>% filter(Leg != "5")

#### FIGURE OF GAM RELATIONSHIPS
inTrain <- createDataPartition(y=df$fgC_cell,
                               times = 1,
                               list = FALSE,
                               p = .8)

training <- df[inTrain,]

gm1 <-gam(fgC_cell ~ s(Temp_C) + Station, data = training)
summary(gm1)

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


# PLOT MODEL
plot_list <- list()

tt = visreg(gm1, xvar = "Temp_C", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("Temp_C" = tt$fit$Temp_C, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
fig1a_line <- line_dat %>% rename(fg_C_cell_mean = "Value", fg_C_cell_95lowerbound = "Lower", fg_C_cell_95upperbound = "Upper") 

dot_dat <- data.frame("Temp_C" = tt$res$Temp_C, "Value" = tt$res$visregRes)
fig1a_points <- dot_dat %>% rename(fg_C_cell = "Value")

x <- densCols(dot_dat$Temp_C,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[1]] <- ggplot() + geom_point(data=dot_dat, aes(x=Temp_C, y = Value, colour = dens), size = 1)+ 
  geom_line(data=line_dat, aes(x=Temp_C, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = Temp_C, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold("Temperature, °C" )))+
  ylab(expression(bold("fg C cell"^{-1}))) + 
  theme_opts

tt = visreg(gm1, xvar = "Station", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("Station" = tt$fit$Station, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
fig1b_line <- line_dat %>% rename(fg_C_cell_mean = "Value", fg_C_cell_95lowerbound = "Lower", fg_C_cell_95upperbound = "Upper") 

dot_dat <- data.frame("Station" = tt$res$Station, "Value" = tt$res$visregRes)
fig1b_points <- dot_dat %>% rename(fg_C_cell = "Value")

x <- densCols(dot_dat$Station ,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[2]] <- ggplot() + geom_point(data=dot_dat, aes(x=Station, y = Value, colour = dens), size = 1)+ 
  geom_point(data=line_dat, aes(x=Station, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_pointrange(data = line_dat, aes(x=Station, y = Value, ymin=Lower, ymax=Upper), col = "red")+
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  scale_x_discrete("Station", breaks = c("001", "070", "140"), labels = c("1", "70", "140"))+
  ylab(expression(bold("fg C cell"^{-1})))  +
  theme_opts

tagss <-c("a", "b")
figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1 , ncol = 2, nrow = 1, legend = "right", common.legend = TRUE)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS8.png", sep = ""), width = 1000, height = 400)


### Save Figure S17 source data
wb <- createWorkbook(type = "xls")
sh1 <- createSheet(wb,sheetName = "Figure S17a_Line")
addDataFrame(data.frame("Figure S17a_Line"=c("Figure S17a_Line")),sheet = sh1,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1a_line,sheet = sh1,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S17a_Points")
addDataFrame(data.frame("Figure S17a_Points"=c("Figure S17a_Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1a_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S17b_Line")
addDataFrame(data.frame("Figure S17b_Line"=c("Figure S17b_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1b_line,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S17b_Points")
addDataFrame(data.frame("Figure S17b Points"=c("Figure S17b Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1b_points,sheet = sh2,row.names = FALSE,startRow = 2)

saveWorkbook(wb,"./source_data/FigureS17.xls")

