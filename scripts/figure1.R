## FIGURE 1: Data distribution and stat model partial residuals 
##          for bacteria abundance

library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggpubr)
library(ggthemes)
library(sp)
library(sf)
library(visreg)
library(caret)
library(colorspace)
library(cowplot)
library(xlsx)

plot_list <- list()

############################################
## Global map of data sites
############################################
curr_bac <- read.csv('./data/curr_bac_final.csv') # Import bacteria abundance data
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
                          plot.margin = unit(c(0,0,0.2,0), "cm")))


r1_gg <- ggplot() + 
  geom_point(data = dat_locs_frame, aes(x=Long, y=Lat), col = "red", size = 0.01)+   geom_polygon(data = land_map, 
                                                      aes(x = X, y = Y, group = group), 
                                                      colour = "black", 
                                                      fill = "black", 
                                                      linewidth = 0.25)+
  theme_opts1 + coord_equal()+ geom_sf(data = Bndry, colour = "black", size = 0.7, fill = NA)


############################################
## Partial residuals plots
############################################

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", linewidth=0.7),
                         axis.text = element_text(size = 14, face = "plain"),
                         axis.title = element_text(size = 16, face = "plain"),
                         plot.margin = unit(c(0.6,0.3,0.3,0.3), "cm")))

#### FIGURE OF FINAL PARAMETRIC MODEL
df <-read.csv('./data/curr_bac_final.csv', header=TRUE)
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
                         panel.border = element_rect(colour = "black", linewidth=0.7),
                         axis.text = element_text(size = 7, face = "plain", colour = "black"),
                         axis.title = element_text(size = 7, face = "plain", colour = "black"),
                         plot.title = element_text(size = 7, face = "plain", colour = "black"),
                         plot.margin = unit(c(0.15,0.1,0.05,0.1), "cm"),
                         legend.key.width = unit(1, "cm"),
                         legend.key.height = unit(2, "cm"),
                         legend.text = element_text(size = 6, face = "plain"),
                         legend.title = element_text(size = 6, face = "plain"),
                         legend.position = 'bottom'))

plot_list <- list()

# PLOT MODEL
tt = visreg(gm1, xvar = "logdepth", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("logdepth" = tt$fit$logdepth, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
fig1b_line <- line_dat %>% rename(log10_prokaryote_abundance_mL_mean = "Value", log10_prokaryote_abundance_mL_95lowerbound = "Lower", log10_prokaryote_abundance_mL_95upperbound = "Upper") 

dot_dat <- data.frame("logdepth" = tt$res$logdepth, "Value" = tt$res$visregRes)
fig1b_points <- dot_dat %>% rename(log10_prokaryote_abundance_mL = "Value") 

x <- densCols(dot_dat$logdepth,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[1]] <- ggplot() + geom_point(data=dot_dat, aes(x=logdepth, y = Value, colour = dens), size = 0.08)+ 
  geom_line(data=line_dat, aes(x=logdepth, y = Value), col = "red", size = 1) + 
  geom_ribbon(data = line_dat, aes(x = logdepth, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette = "BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  theme_bw()+
  scale_x_continuous(name = expression(paste("Depth, m", sep = "")), breaks = c(0,1,2,3), labels = c(1,10,100,1000)) + 
  ylab(expression(paste("log"[10], "(Abundance, mL"^-1, ")", sep = "")))+
  theme_opts + theme(legend.position = "none")

tt = visreg(gm1, xvar = "log_nitrate", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("log_nitrate" = tt$fit$log_nitrate, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
fig1c_line <- line_dat %>% rename(log10_prokaryote_abundance_mL_mean = "Value", log10_prokaryote_abundance_mL_95lowerbound = "Lower", log10_prokaryote_abundance_mL_95upperbound = "Upper") 

dot_dat <- data.frame("log_nitrate" = tt$res$log_nitrate, "Value" = tt$res$visregRes)
fig1c_points <- dot_dat %>% rename(log10_prokaryote_abundance_mL = "Value")

x <- densCols(dot_dat$log_nitrate,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[2]] <- ggplot() + geom_point(data=dot_dat, aes(x=log_nitrate, y = Value, colour = dens), size = 0.08)+ 
  geom_line(data=line_dat, aes(x=log_nitrate, y = Value), col = "red", size = 1) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = log_nitrate, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette = "BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  scale_x_continuous(name = expression(paste("Nitrate ", mu, "mol kg"^-1,sep = "")), breaks = c(-2,-1,0,1), labels = c(0.01, 0.1, 1, 10)) + 
  ylab(expression(paste("log"[10], "(Abundance, mL"^-1, ")", sep = ""))) +
  ylim(c(4.5,6.5)) + 
  theme_opts+ theme(legend.position = "none")

tt = visreg(gm1, xvar = "aou", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("aou" = tt$fit$aou, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
fig1d_line <- line_dat %>% rename(log10_prokaryote_abundance_mL_mean = "Value", log10_prokaryote_abundance_mL_95lowerbound = "Lower", log10_prokaryote_abundance_mL_95upperbound = "Upper") 

dot_dat <- data.frame("aou" = tt$res$aou, "Value" = tt$res$visregRes)
fig1d_points <- dot_dat %>% rename(log10_prokaryote_abundance_mL = "Value")

x <- densCols(dot_dat$aou,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[3]] <- ggplot() + geom_point(data=dot_dat, aes(x=aou, y = Value, colour = dens), size = 0.08)+ 
  geom_line(data=line_dat, aes(x=aou, y = Value), col = "red", size = 1) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = aou, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette = "BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(paste("AOU, ", mu, "mol kg"^-1, sep = ""))) + 
  ylab(expression(paste("log"[10], "(Abundance, mL"^-1, ")", sep = ""))) +
  ylim(c(4.5,6.5))+ 
  theme_opts+theme(legend.position = "none")

tt = visreg(gm1, xvar = "temperature", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("temperature" = tt$fit$temperature, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
fig1e_line <- line_dat %>% rename(log10_prokaryote_abundance_mL_mean = "Value", log10_prokaryote_abundance_mL_95lowerbound = "Lower", log10_prokaryote_abundance_mL_95upperbound = "Upper") 

dot_dat <- data.frame("temperature" = tt$res$temperature, "Value" = tt$res$visregRes)
fig1e_points <- dot_dat %>% rename(log10_prokaryote_abundance_mL = "Value")

x <- densCols(dot_dat$temperature ,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[4]] <- ggplot() + geom_point(data=dot_dat, aes(x=temperature, y = Value, colour = dens), size = 0.08)+ 
  geom_line(data=line_dat, aes(x=temperature, y = Value), col = "red", size = 1) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = temperature, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette = "BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression("Temperature, °C" )) + 
  ylab(expression(paste("log"[10], "(Abundance, mL"^-1, ")", sep = ""))) +
  ylim(c(4.5,6.5))+ 
  theme_opts+theme(legend.position = "none")

tt = visreg(gm1, xvar = "logchlo", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("logchlo" = tt$fit$logchlo, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
fig1f_line <- line_dat %>% rename(log10_prokaryote_abundance_mL_mean = "Value", log10_prokaryote_abundance_mL_95lowerbound = "Lower", log10_prokaryote_abundance_mL_95upperbound = "Upper") 

dot_dat <- data.frame("logchlo" = tt$res$logchlo, "Value" = tt$res$visregRes)
fig1f_points <- dot_dat %>% rename(log10_prokaryote_abundance_mL = "Value")

x <- densCols(dot_dat$logchlo ,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[5]] <- ggplot() + geom_point(data=dot_dat, aes(x=logchlo, y = Value, colour = dens), size = 0.08)+ 
  geom_line(data=line_dat, aes(x=logchlo, y = Value), col = "red", size = 1) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = logchlo, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette ="BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  scale_x_continuous(name = expression(paste("Chlorophyll, mg m"^-3, sep = "")), breaks = c(-1,0,1), labels = c(0.1,1,10)) + 
  ylab(expression(paste("log"[10], "(Abundance, mL"^-1, ")", sep = ""))) +
  ylim(c(4.5,6.5))+
  theme_opts + theme(legend.position = "none")

plot_list[[6]] <- ggplot(data=dot_dat, aes(x=logchlo, y = Value, colour = dens))+
  geom_point()+
  lims(x = c(0,0), y = c(0,0))+
  theme_void()+
  theme(legend.position = c(0.5,0.6),
        legend.key.size = unit(0.25, "cm"),
        legend.key.height = unit(0.35, "cm"),
        legend.text = element_text(size = 7, face = "plain"),
        legend.title = element_text(size = 7, face = "plain"),
        legend.key.width = unit(0.3, "cm"))+
  scale_colour_continuous_sequential(palette ="BluYl", rev = FALSE, name = "Number of\nsamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  ">240"))

kk = plot_grid(r1_gg, 
          plot_grid(plotlist = plot_list, ncol = 2, labels = c("b", "c", "d", "e", "f"), label_y = 1.05, label_size = 7),
          ncol = 1, rel_heights = c(0.33,1), rel_widths = c(2.5,1), labels = c('a', ' '), label_size = 7)

ggsave(kk, filename = paste("./figures/Figure1.pdf", sep = ""), width = 88, height = 150, units = "mm")



### Save Figure 1 source data
dat_locs <- dat_locs %>% dplyr::rename(Longitude = "long", Latitude = "lat")
wb <- createWorkbook(type = "xls")
sh1 <- createSheet(wb,sheetName = "Figure 1a")
addDataFrame(data.frame("Figure 1a"=c("Figure 1a")),sheet = sh1,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(dat_locs,sheet = sh1,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure 1b_Points")
addDataFrame(data.frame("Figure 1b_Points"=c("Figure 1b_Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1b_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure 1b_Line")
addDataFrame(data.frame("Figure 1b_Line"=c("Figure 1b_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1b_line,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure 1c_Points")
addDataFrame(data.frame("Figure 1c"=c("Figure 1c")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1c_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure 1c_Line")
addDataFrame(data.frame("Figure 1c_Line"=c("Figure 1c_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1c_line,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure 1d_Points")
addDataFrame(data.frame("Figure 1d_Points"=c("Figure 1d_Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1d_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure 1d_Line")
addDataFrame(data.frame("Figure 1d_Line"=c("Figure 1d_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1d_line,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure 1e_Points")
addDataFrame(data.frame("Figure 1e_Points"=c("Figure 1e_Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1e_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure 1e_Line")
addDataFrame(data.frame("Figure 1e_Line"=c("Figure 1e_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1e_line,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure 1f_Points")
addDataFrame(data.frame("Figure 1f_Points"=c("Figure 1f_Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1f_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure 1f_Line")
addDataFrame(data.frame("Figure 1f_Line"=c("Figure 1f_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(fig1f_line,sheet = sh2,row.names = FALSE,startRow = 2)

saveWorkbook(wb,"./source_data/Figure1/Figure1.xls")


