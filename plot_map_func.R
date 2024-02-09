library(raster)
library(maptools)
library(ggplot2)
library(scales)
library(rgdal)
library(sf)
library(colorspace)

map_figure <- function(rel_map, curr_pal = "YlOrBr", legend_name, plot_lims = NA, unlog_leg = FALSE, dp = 1, var_name = NA){
  ## Import land values for ggplot
  land_map <- read.csv('./scripts/land_shapefiles/robinson_proj.csv')

  rel_map_curr <- matrix(unlist(as.vector(rel_map[,3])), nrow = 360, ncol = 180)
  
  rel_map <- rel_map_curr[c(181:360,1:180),180:1]
  # Set up the raster
  e <- extent(c(-180, 180, -90, 90))
  
  r2 <- raster(nrows = 180, ncols = 360, ext = e)
  r2[] <- t(rel_map)
  
  # Set up projections   
  robinson_proj <- CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  cc_rob <-projectRaster(r2, crs = robinson_proj, over = T)
  
  # To convert your RasterLayer to a data.frame, you need to convert it to
  # a SpatialPixelsDataFrame first
  r.1 <- rasterToPoints(cc_rob)
  r.1.df <- as.data.frame(r.1)
  colnames(r.1.df) <-c("x","y","Sp")

  
  num_dp = dp # Number of dp to round legend key
  
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
                            plot.background = element_rect(fill="white"),
                            panel.border = element_blank(),
                            axis.line = element_blank(),
                            axis.text.x = element_blank(),
                            axis.text.y = element_blank(),
                            axis.ticks = element_blank(),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            plot.tag = element_text(face = "bold", size = 12),
                            plot.tag.position = "topleft",
                            plot.title = element_text(size = 7, hjust = 0.5, face = "bold"),
                            plot.subtitle=element_blank(),
                            legend.text = element_text(size = 7),
                            legend.title = element_text(size = 6),
                      legend.key.width = unit(0.4, "cm"),
                      legend.key.height = unit(0.2, "cm"),
                      legend.box.spacing = unit(0, "lines"),
                      legend.position = "bottom",
                      legend.box.margin = unit(c(0,0,0,0), "lines"),
                      plot.margin = unit(c(0,0,0,0), "lines")))
  
  
  if(is.na(plot_lims[1]) == TRUE){
    plot_limits = c(min(r.1.df$Sp, na.rm=TRUE), max(r.1.df$Sp, na.rm=TRUE))
  }else{plot_limits = plot_lims}
  
  r1_gg <- ggplot(r.1.df) + 
    geom_tile(aes(x=x, y=y, fill = Sp))+   geom_polygon(data = land_map, 
                                              aes(x = X, y = Y, group = group), 
                                              colour = "black", 
                                              fill = "black", 
                                              linewidth = 0.1)+
    scale_fill_continuous_sequential(palette = curr_pal, rev = TRUE, limits = c(plot_limits[1], plot_limits[2]),oob = squish, guide = guide_colorbar(title = legend_name),
                                     name = "", na.value='lightgrey',breaks = seq(plot_limits[1],plot_limits[2],length.out = 3),
                         labels = format(round(seq(plot_limits[1],plot_limits[2],length.out = 3),digits = num_dp), nsmall = num_dp))  +
    theme_opts1  + coord_equal(ratio = 1.2)+ ggtitle(var_name) + geom_sf(data = Bndry, colour = "black", linewidth = 0.5, fill = NA)
  
  if(unlog_leg == TRUE){
    r1_gg <- r1_gg +     scale_fill_continuous_sequential(palette = curr_pal, rev = TRUE, limits = c(plot_limits[1], plot_limits[2]),oob = squish, guide = guide_colorbar(title = legend_name),
                                                          name = "", na.value='lightgrey',breaks = seq(plot_limits[1],plot_limits[2],length.out = 3),
                                                          labels = format(round(10^seq(plot_limits[1],plot_limits[2],length.out = 3),digits = num_dp), nsmall = num_dp))
  }
  
  return(r1_gg)
}



map_figure_delta <- function(rel_map, curr_pal = "Red-Blue", legend_name, plot_lims = NA, var_name = NA){
  ## Import land values for ggplot
  land_map <- read.csv('./scripts/land_shapefiles/robinson_proj.csv')
  
  rel_map <- rel_map[c(181:360,1:180),180:1]
  # Set up the raster
  e <- extent(c(-180, 180, -90, 90))
  
  r2 <- raster(nrows = 180, ncols = 360, ext = e)
  r2[] <- t(rel_map)
  
  # Set up projections   
  robinson_proj <- CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  cc_rob <-projectRaster(r2, crs = robinson_proj, over = T)
  
  # To convert your RasterLayer to a data.frame, you need to convert it to
  # a SpatialPixelsDataFrame first
  r.1 <- rasterToPoints(cc_rob)
  r.1.df <- as.data.frame(r.1)
  colnames(r.1.df) <-c("x","y","Sp")
  
  
  num_dp = 2 # Number of dp to round legend key
  
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
                            plot.background = element_rect(fill="white"),
                            panel.border = element_blank(),
                            plot.margin = unit(c(0,0,0,0),"lines"),
                            axis.line = element_blank(),
                            axis.text.x = element_blank(),
                            axis.text.y = element_blank(),
                            axis.ticks = element_blank(),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            plot.tag = element_text(face = "bold", size = 18),
                            plot.tag.position = "topleft",
                            plot.title = element_text(size = 20, face = "bold",hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                            plot.subtitle=element_text(size= 12, color="black"),
                            legend.text = element_text(size = 12, face = "bold")))
  
  
  if(is.na(plot_lims[1]) == TRUE){
    plot_limits = c(min(r.1.df$Sp, na.rm=TRUE), max(r.1.df$Sp, na.rm=TRUE))
  }else{plot_limits = plot_lims}
  
  r1_gg <- ggplot(r.1.df, aes(x=x, y=y)) + 
    geom_tile(aes(fill = Sp))+   geom_polygon(data = land_map, 
                                              aes(x = X, y = Y, group = group), 
                                              colour = "black", 
                                              fill = "black", 
                                              size = 0.25)+
    # scale_fill_gradient2(low = 'royalblue1', mid = "white",  high = 'red', limits = c(plot_limits[1], plot_limits[2]),oob = squish, guide = guide_colorbar(label.position = "left", title.hjust = 1, title = expression(paste('kg C m' ^ '-2'))),
    scale_fill_continuous_diverging(palette = curr_pal, rev = TRUE, limits = c(plot_limits[1], plot_limits[2]),oob = squish, guide = guide_colorbar(title = legend_name, title.theme = element_text(face = "bold")),
                                     name = "", na.value='lightgrey',breaks = seq(plot_limits[1],plot_limits[2],length.out = 5),
                                     labels = format(round(seq(plot_limits[1],plot_limits[2],length.out = 5),digits = num_dp), nsmall = num_dp))  +
    theme_opts1  + coord_equal()+ ggtitle(var_name) 
  #+ geom_sf(data = Bndry, colour = "black", size = 0.7, fill = NA)
  
  return(r1_gg)
}

