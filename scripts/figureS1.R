## SUPPLEMENTARY FIGURE S1: Data distribution and stat model residuals for bacteria cell carbon model

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggthemes)
library(sp)
library(viridis)
library(lme4)
library(colorspace)
library(cowplot)
library(caret)
library(visreg)
library(xlsx)

plot_list <- list()

############################################
## Global map of data sites
############################################
curr_bac <- read.csv('./data/malaspina_locations.csv') # Import malaspina locations
dat_locs <- curr_bac %>% dplyr::select(Lon, Lat)  %>% na.omit
dat_locs[dat_locs$Lon < 0,"Lon"] <- dat_locs[dat_locs$Lon < 0,"Lon"] + 360

## Import land values for ggplot
land_map <- read.csv('./scripts/land_shapefiles/robinson_proj.csv')
robinson_proj <- CRS("+proj=robin +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

dat_locs_sp <- SpatialPoints(dat_locs)
crs(dat_locs_sp) <- "+proj=longlat +datum=WGS84"

dat_locs_rob <- spTransform(dat_locs_sp, robinson_proj)

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
  geom_point(data = dat_locs_frame, aes(x=Long, y=Lat), col = "red", size = 1.5)+   
  geom_polygon(data = land_map,  aes(x = X, y = Y, group = group),  colour = "black", 
               fill = "black", linewidth = 0.25)+
  theme_opts1 + coord_equal()+ geom_sf(data = Bndry, colour = "black", size = 0.7, fill = NA)

############################################
## Partial residuals plots
############################################

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", size=0.7),
                         axis.text = element_text(size = 14, face = "bold"),
                         axis.title = element_text(size = 16, face = "bold"),
                         plot.margin = unit(c(0.6,0.3,0.3,0.3), "cm")))


#### FIGURE OF FINAL PARAMETRIC MODEL
df <-read.csv('./data/raw_observations/prokaryote_cellcarbon.csv', header=TRUE)
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

#Split the data set, sample, into a training data set (80%) and test data set (20%).
#The training data is used to train the model, and then the model is run on the test data to check how well it is classifying the data.

inTrain <- createDataPartition(y=df$fgC_cell,
                               times = 1,
                               list = FALSE,
                               p = .8)

training <- df[inTrain,]
testing <- df[-inTrain,]

### GLMM WITH TEMPERATURE AND RANDOM EFFECT FOR STATION
gm1 <- lmer(fgC_cell ~ poly(Temp_C,3)+(1|Station), data = training)
summary(gm1)

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", size=0.7),
                         axis.text = element_text(size = 16, face = "bold", colour = "black"),
                         axis.title = element_text(size = 18, face = "bold", colour = "black"),
                         plot.title = element_text(size = 16, face = "bold.italic", colour = "black"),
                         plot.margin = unit(c(0.6,0.3,0.3,0.3), "cm"),
                         legend.key.width = unit(1, "cm"),
                         legend.key.height = unit(1., "cm"),
                         legend.text = element_text(size = 20, face = "bold"),
                         legend.title = element_text(size = 20, face = "bold"),
                         legend.position = 'right'))

plot_list <- list()
# PLOT MODEL
tt = visreg(gm1, xvar = "Temp_C", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("Temp_C" = tt$fit$Temp_C, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
figS1b_line <- line_dat %>% rename(fg_C_cell = "Value", fg_C_cell_95lowerbound = "Lower", fg_C_cell_95upperbound = "Upper") 

dot_dat <- data.frame("Temp_C" = tt$res$Temp_C, "Value" = tt$res$visregRes)
figS1b_points <- dot_dat %>% rename(fg_C_cell = "Value")

x <- densCols(dot_dat$Temp_C,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[1]] <- ggplot() + geom_point(data=dot_dat, aes(x=Temp_C, y = Value, colour = dens), size = 1)+ 
  geom_line(data=line_dat, aes(x=Temp_C, y = Value), col = "red", size = 2) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = Temp_C, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_continuous_sequential(palette = "BluYl", rev = FALSE, name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold("Temperature, °C" )))+
  ylab(expression(bold("fg C cell"^{-1}))) + 
  theme_opts

kk = plot_grid(r1_gg, plot_list[[1]], ncol = 2, labels = c("a", "b"))
ggpubr::ggexport(kk, filename ="./figures/FigureS1.png",  width = 800, height = 250)

### Save Figure S1 source data
dat_locs <- dat_locs %>% dplyr::rename(Longitude = "Lon", Latitude = "Lat")
wb <- createWorkbook(type = "xls")
sh1 <- createSheet(wb,sheetName = "Figure S1a")
addDataFrame(data.frame("Figure S1a"=c("Figure S1a")),sheet = sh1,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(dat_locs,sheet = sh1,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S1b_Points")
addDataFrame(data.frame("Figure S1b_Points"=c("Figure S1b_Points")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(figS1b_points,sheet = sh2,row.names = FALSE,startRow = 2)

sh2 <- createSheet(wb,sheetName = "Figure S1b_Line")
addDataFrame(data.frame("Figure S1b_Line"=c("Figure S1b_Line")),sheet = sh2,row.names = FALSE,startRow = 1,col.names = FALSE)
addDataFrame(figS1b_line,sheet = sh2,row.names = FALSE,startRow = 2)

saveWorkbook(wb,"./source_data/FigureS1.xls")
