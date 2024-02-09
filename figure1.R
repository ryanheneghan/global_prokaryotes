## FIGURE 1: Data distribution and stat model residuals for bacteria abundance

rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggpubr)
library(ggthemes)
library(sp)
library(rgdal)
library(randomForest)
library(visreg)
library(caret)
library(colorspace)
library(cowplot)

setwd("~/Desktop/Papers/Bacteria_Census/")

source('./scripts/plot_map_func.R')

plot_list <- list()

############################################
## Global map of data sites
############################################
curr_bac <- read.csv('./curr_bac_final.csv') # Import bacteria abundance data
dat_locs <- curr_bac %>% dplyr::select(long, lat)  %>% na.omit
dat_locs[dat_locs$long < 0,"long"] <- dat_locs[dat_locs$long < 0,"long"] + 360

## Import land values for ggplot
land_map <- read.csv('./scripts/land_shapefiles/robinson_proj.csv')
robinson_proj <- CRS("+proj=robin +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

dat_locs_sp <- SpatialPoints(dat_locs)
crs(dat_locs_sp) <- "+proj=longlat +datum=WGS84"

dat_locs_rob <- spTransform(dat_locs_sp, robinson_proj)
plot(dat_locs_rob)

dat_locs_frame <- as.data.frame(dat_locs_rob)
colnames(dat_locs_frame) <- c("Long", "Lat")

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


r1_gg <- ggplot() + 
  geom_point(data = dat_locs_frame, aes(x=Long, y=Lat), col = "red", size = 1)+   geom_polygon(data = land_map, 
                                                      aes(x = X, y = Y, group = group), 
                                                      colour = "black", 
                                                      fill = "black", 
                                                      size = 0.25)+
  theme_opts1 + coord_equal()+ geom_sf(data = Bndry, colour = "black", size = 0.7, fill = NA)


############################################
## Partial residuals plots
############################################

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", linewidth=0.7),
                         axis.text = element_text(size = 14, face = "bold"),
                         axis.title = element_text(size = 16, face = "bold"),
                         plot.margin = unit(c(0.6,0.3,0.3,0.3), "cm")))

#### FIGURE OF FINAL PARAMETRIC MODEL
df <-read.csv('./curr_bac_final.csv', header=TRUE)
df <- df %>% dplyr::filter(silicate > 0 & phosphate > 0 & nitrate > 0 & oxygen > 0)

df$log_phosphate <- log10(df$phosphate)
df$log_nitrate <- log10(df$nitrate)
df$log_silicate <- log10(df$silicate)

df1 <- df %>% #dplyr::filter(10^logdepth > z_eu) %>% 
  dplyr::mutate(si_star = silicate - nitrate,
                n_star = nitrate-16*phosphate) %>%
  dplyr::select(c("logabund","logchlo", "logdepth","temperature", "oxygen", "log_silicate", "log_nitrate", "log_phosphate", "n_star", "si_star", "aou"))%>%
  dplyr::filter(n_star < 50)

inTrain <- createDataPartition(y=df1$logabund,
                               times = 1,
                               list = FALSE,
                               p = .8)

df1 <- df1 [inTrain,]

gm1 <- lm(logabund ~ poly(logdepth,5) + log_nitrate + poly(aou, 2) + temperature + logchlo, data = df1)
summary(gm1)



theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", size=0.7),
                         axis.text = element_text(size = 16, face = "bold", colour = "black"),
                         axis.title = element_text(size = 18, face = "bold", colour = "black"),
                         plot.title = element_text(size = 16, face = "bold.italic", colour = "black"),
                         plot.margin = unit(c(0.6,0.3,0.3,0.3), "cm"),
                         legend.key.width = unit(1, "cm"),
                         legend.key.height = unit(2, "cm"),
                         legend.text = element_text(size = 20, face = "bold"),
                         legend.title = element_text(size = 20, face = "bold"),
                         legend.position = 'bottom'))

plot_list <- list()
# PLOT MODEL
tt = visreg(gm1, xvar = "logdepth", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("logdepth" = tt$fit$logdepth, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("logdepth" = tt$res$logdepth, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$logdepth,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[1]] <- ggplot() + geom_point(data=dot_dat, aes(x=logdepth, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=logdepth, y = Value), col = "red", size = 1.5) + 
  geom_ribbon(data = line_dat, aes(x = logdepth, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette = "BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  theme_bw()+
  scale_x_continuous(name = expression(bold(paste("Depth, m", sep = ""))), breaks = c(0,1,2,3), labels = c(1,10,100,1000)) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = ""))))+
  theme_opts + theme(legend.position = "none")

tt = visreg(gm1, xvar = "log_nitrate", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("log_nitrate" = tt$fit$log_nitrate, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("log_nitrate" = tt$res$log_nitrate, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$log_nitrate,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[2]] <- ggplot() + geom_point(data=dot_dat, aes(x=log_nitrate, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=log_nitrate, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = log_nitrate, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette = "BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  scale_x_continuous(name = expression(bold(paste("Nitrate ", mu, "mol kg"^-1,sep = ""))), breaks = c(-3,-2,-1,0,1), labels = c(0.001, 0.01, 0.1, 1, 10)) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5)) + 
  theme_opts+ theme(legend.position = "none")

tt = visreg(gm1, xvar = "aou", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("aou" = tt$fit$aou, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("aou" = tt$res$aou, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$aou,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[3]] <- ggplot() + geom_point(data=dot_dat, aes(x=aou, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=aou, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = aou, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette = "BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold(paste("AOU (", mu, "mol kg"^-1, ")", sep = "")))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+ 
  theme_opts+theme(legend.position = "none")

tt = visreg(gm1, xvar = "temperature", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("temperature" = tt$fit$temperature, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("temperature" = tt$res$temperature, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$temperature ,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[4]] <- ggplot() + geom_point(data=dot_dat, aes(x=temperature, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=temperature, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = temperature, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette = "BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold("Temperature, °C" ))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+ 
  theme_opts+theme(legend.position = "none")

tt = visreg(gm1, xvar = "logchlo", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("logchlo" = tt$fit$logchlo, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("logchlo" = tt$res$logchlo, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$logchlo ,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[5]] <- ggplot() + geom_point(data=dot_dat, aes(x=logchlo, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=logchlo, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = logchlo, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette ="BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  scale_x_continuous(name = expression(bold(paste("Chlorophyll, mg m"^-3, sep = ""))), breaks = c(-1,0,1), labels = c(0.1,1,10)) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+
  theme_opts + theme(legend.position = "none")

plot_list[[6]] <- ggplot(data=dot_dat, aes(x=logchlo, y = Value, colour = dens))+
  geom_point()+
  lims(x = c(0,0), y = c(0,0))+
  theme_void()+
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(1, "cm"),
        legend.key.height = unit(1.2, "cm"),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.key.width = unit(1, "cm"))+
  scale_colour_continuous_sequential(palette ="BluYl", rev = FALSE, name = "Number of\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))

kk = plot_grid(r1_gg, 
          plot_grid(plotlist = plot_list, ncol = 2, labels = c("b", "c", "d", "e", "f"), label_size = 20),
          ncol = 1, rel_heights = c(0.33,1), rel_widths = c(2.5,1), labels = c('a', ' '), label_size = 20)

ggpubr::ggexport(kk, filename = paste("./figures/Figure1.png", sep = ""), width = 700, height = 1200)

