rm(list = ls())

library(dplyr)
library(tidyr)
library(viridis)
library(ggpointdensity)
library(ggplot2)
library(ggpubr)
library(lme4)
library(visreg)
library(mgcv)

setwd("~/Desktop/Papers/Bacteria_Census/growth_data")

source('~/Desktop/Papers/Bacteria_Census/scripts/plot_map_func.R')

df <- read.csv("all_sgr.csv")
df <- df %>% dplyr::mutate(log10SGR = log10(SGR),
                           log10CHL = log10(CHL),
                           log10Depth = log10(Depth))

gge_dat <- read.csv("robinson_respiration_clean.csv")
############################################
## EXTENDED DATA FIGURE 3 global maps of SGR and GGE data sites
dat_locs <- df %>% dplyr::select(Longitude, Latitude)  %>% na.omit
dat_locs_gge <- gge_dat %>% dplyr::select(Long, Lat) %>% na.omit

## Import land values for ggplot
land_map <- read.csv('~/Desktop/Papers/Bacteria_Census/scripts/land_shapefiles/robinson_proj.csv')
robinson_proj <- CRS("+proj=robin +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

dat_locs_sp <- SpatialPoints(dat_locs)
dat_locs_sp2 <- SpatialPoints(dat_locs_gge)

crs(dat_locs_sp) <- "+proj=longlat +datum=WGS84"
crs(dat_locs_sp2) <- "+proj=longlat +datum=WGS84"

dat_locs_rob <- spTransform(dat_locs_sp, robinson_proj)
dat_locs_rob2 <- spTransform(dat_locs_sp2, robinson_proj)

dat_locs_frame <- as.data.frame(dat_locs_rob)
colnames(dat_locs_frame) <- c("Long", "Lat")

dat_locs_frame2 <- as.data.frame(dat_locs_rob2)
colnames(dat_locs_frame2) <- c("Long", "Lat")

## Boundary layer for map
Bndry <- tibble(x = seq(-180, 180, by = 1), y = -90) %>%
  bind_rows(tibble(x = 180, y = seq(-90, 90, by = 1))) %>%
  bind_rows(tibble(x = seq(180, -180, by = -1), y = 90)) %>%
  bind_rows(tibble(x = -180, y = seq(90, -90, by = -1))) %>%
  as.matrix() %>%
  list() %>%
  st_polygon() %>%
  st_sfc(crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform(crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# 
theme_opts1 <- list(theme(panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.background = element_blank(),
                          panel.border = element_blank(),
                          axis.line = element_blank(),
                          axis.text.x = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks = element_blank(),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          plot.title = element_text(size=12, hjust = 0.5),
                          legend.key.width = unit(2, "cm"),
                          plot.margin = unit(c(-0.2,-0.5,-0.5,-0.5), "cm")))

plot_list <- list()
plot_list[[1]] <- ggplot() + 
  geom_point(data = dat_locs_frame, aes(x=Long, y=Lat), col = "red", size = 1)+   geom_polygon(data = land_map, 
                                                                                               aes(x = X, y = Y, group = group), 
                                                                                               colour = "black", 
                                                                                               fill = "black", 
                                                                                               size = 0.25)+
  theme_opts1 + coord_equal()+ geom_sf(data = Bndry, colour = "black", size = 0.7, fill = NA)


plot_list[[2]] <- ggplot() + 
  geom_point(data = dat_locs_frame2, aes(x=Long, y=Lat), col = "red", size = 1)+   geom_polygon(data = land_map, 
                                                                                               aes(x = X, y = Y, group = group), 
                                                                                               colour = "black", 
                                                                                               fill = "black", 
                                                                                               size = 0.25)+
  theme_opts1 + coord_equal()+ geom_sf(data = Bndry, colour = "black", size = 0.7, fill = NA)

tagss <-c("a", "b")
figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1.5, ncol = 1, nrow =2, legend = "right", common.legend = TRUE)
ggpubr::ggexport(figure, filename = paste("~/Desktop/Papers/Bacteria_Census/figures/FigureS10.png", sep = ""), width = 400, height =450)


#### FIGURE OF FINAL PARAMETRIC MODEL
all_data <- read.csv("all_sgr.csv")
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
# PLOT MAIN EFFECTS
plot_list <- list()
plot_list <- list()

tt = visreg(lmer_temp, xvar = "SST", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("Temp_C" = tt$fit$SST, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("Temp_C" = tt$res$SST, "Value" = tt$res$visregRes)

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
dot_dat <- data.frame("log10CHL" = tt$res$log10CHL, "Value" = tt$res$visregRes)

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
ggpubr::ggexport(figure, filename = paste("~/Desktop/Papers/Bacteria_Census/figures/FigureS13.png", sep = ""), width = 1000, height = 400)


#### REGIONAL MODEL ASSESSMENT
df <-read.csv('./curr_bac_final.csv', header=TRUE)
df <- df %>% dplyr::filter(silicate > 0 & phosphate > 0 & nitrate > 0 & oxygen > 0)
# From Morel et al. 2007 Examining the consistency of products derived from various ocean colorsensors in open ocean (Case 1) waters in the perspective of a multi-sensor approach
#df <- df %>% dplyr::mutate(z_eu = 10^(1.524-0.460*logchlo-0.00051*logchlo^2+0.0282*logchlo^3))

df$log_phosphate <- log10(df$phosphate)
df$log_nitrate <- log10(df$nitrate)
df$log_silicate <- log10(df$silicate)

df1 <- df %>% #dplyr::filter(10^logdepth > z_eu) %>% 
  dplyr::mutate(si_star = silicate - nitrate,
                n_star = nitrate-16*phosphate) %>%
  dplyr::select(c("lat","long","logabund","logchlo", "logdepth","temperature", "oxygen", "log_silicate", "log_nitrate", "log_phosphate", "n_star", "si_star", "aou"))%>%
  dplyr::filter(n_star < 50)

inTrain <- createDataPartition(y=df1$logabund,
                               times = 1,
                               list = FALSE,
                               p = .8)
all1 <- df1
training1 <- df1 [inTrain,]

gm1 <- lm(logabund ~ poly(logdepth,5) + temperature + log_nitrate + logchlo + poly(aou,2), data = training1)

min_lat <- c(0,0,30,60)
max_lat <- c(90,30,60,90)
regions <- c("Global", "Polar", "Temperate", "Tropical")

#Function to calculate mean square error
rmse <- function(error)
{
  sqrt(mean(error)^2)
}

#Calculate mean absolute error
mae <- function(error)
{
  mean(abs(error))
}

plot_list <- list()

for(i in 1:length(regions)){
  reg_dat1 <- all1 %>% dplyr::filter(abs(lat) > min_lat[i] & abs(lat) <= max_lat[i])
  x_test1 <- reg_dat1 %>% dplyr::select(logdepth,logchlo, temperature, log_nitrate,aou)
  reg_pred1 <- predict(gm1, x_test1, type = "response")
  error1<- (reg_pred1 - reg_dat1$logabund)
  
  
  gg_dat <- data.frame("x_dat" = reg_dat1$logabund, "y_dat" = reg_pred1)
  curr_r2 <- cor(gg_dat$x_dat, gg_dat$y_dat)^2
  ann_label <- expression("R"^2 == curr_r2)
  
  plot_list[[i]] <- ggplot(gg_dat) + geom_point(aes(x=x_dat,y=y_dat), shape = 1) + geom_abline(slope = 1, intercept = 0, colour = "red", size = 1) +
    theme_bw() + theme_opts + xlab(expression(bold(paste("Observed log"[10], "(cells ml"^-1, ")", sep = "")))) +
    ylab(expression(bold(paste("Predicted log"[10], "(cells ml"^-1, ")", sep = "")))) +
    labs(title = paste(regions[i], ", n = ", length(reg_pred1), sep = "")) + xlim(c(min(df1$logabund),max(df1$logabund))) + ylim(c(4.3,6.2))+
    annotate("label", x =6.5, y = 4.4, parse = TRUE, label =paste("R ^ 2 == ", round(curr_r2*100,1),sep = ""), size = 6, fontface = "bold")
}


tagss <-c("A)", "B)", "C)", "D)")
figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1 , ncol = 2, nrow = 2)
ggpubr::ggexport(figure, filename = paste("./figures/supp_figure4b.png", sep = ""), width = 800, height = 700)

