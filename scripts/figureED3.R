## EXTENDED DATA FIGURE 3 global maps of SGR and GGE data sites

library(dplyr)
library(tidyr)
library(viridis)
library(ggpointdensity)
library(ggplot2)
library(ggpubr)


df <- read.csv("./data/all_sgr.csv")
df <- df %>% dplyr::mutate(log10SGR = log10(SGR),
                           log10CHL = log10(CHL),
                           log10Depth = log10(Depth))

gge_dat <- read.csv("./data/robinson_respiration_clean.csv")

############################################
dat_locs <- df %>% dplyr::select(Longitude, Latitude)  %>% na.omit
dat_locs_gge <- gge_dat %>% dplyr::select(Long, Lat) %>% na.omit

## Import land values for ggplot
land_map <- read.csv('./scripts/land_shapefiles/robinson_proj.csv')
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
ggpubr::ggexport(figure, filename = paste("./figures/FigureED3.png", sep = ""), width = 400, height =450)

