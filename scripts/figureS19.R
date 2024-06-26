## Figure S19: Raw relationships for SPR

library(dplyr)
library(tidyr)
library(viridis)
library(ggpointdensity)
library(ggplot2)
library(ggpubr)


df <- read.csv('./data/raw_observations/prokaryote_sgr.csv', header=TRUE)
df <- df %>% dplyr::mutate(log10SGR = log10(SGR),
                           log10CHL = log10(CHL),
                           log10Depth = log10(Depth))

x_names <- c("SST", "log10CHL","source")
x_save_name <- c(expression(bold("Temperature, °C" )),
                 expression(bold(paste("log"[10], "(Chlorophyll, mg m"^{-3} ,")", sep = ""))),
                 expression(bold("Dataset")))

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

tagss <-c("a", "b", "c")
curr_plot <- list()
wb <- createWorkbook(type = "xls")

for(i in 1:length(x_names)){
  curr_dat <- df[,c("log10SGR",x_names[i])]
  source_data <- curr_dat
  names(curr_dat) <- c("y_var", "x_var")
  x <- densCols(curr_dat$x_var,curr_dat$y_var, colramp=colorRampPalette(c("black", "white")))
  
  ## Use densCols() output to get density at each point
  curr_dat$dens <- col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                              "#FCFF00", "#FF9400", "#FF3100"))(256)
  curr_dat$col <- cols[curr_dat$dens]
  
  ## Reorder rows so that densest points are plotted on top
  curr_dat <- curr_dat[order(curr_dat$dens),]
  
  curr_plot[[i]] <- ggplot(curr_dat) +
    geom_point(aes(x=x_var,y=y_var, colour = dens))+
    theme_bw() + scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
    ylab(expression(bold("fg C cell"^{-1})))+ 
    xlab(x_save_name[i]) + guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) +
    theme_opts
  
  if(x_names[i] == "source"){
    curr_plot[[i]] <- curr_plot[[i]] + scale_x_discrete(x_save_name[i], breaks = c("Blanes_Bay", "HotMix", "Latitud", "Malaspina", "Kirchman_2009"), labels = c("Blanes Bay", "HotMix", "Latitud", "Malaspina-2010", "Polar"))
  }
  
  ## Save source data for Figure S16
  header_name <- paste("Figure S19", tagss[i], sep = "")
  
  sh1 <- createSheet(wb,sheetName = header_name)
  addDataFrame(data.frame(header_name=c(header_name)),sheet = sh1,row.names = FALSE,startRow = 1,col.names = FALSE)
  addDataFrame(source_data,sheet = sh1,row.names = FALSE,startRow = 2)
}


saveWorkbook(wb,"./source_data/FigureS19.xls")

figure <- ggpubr::ggarrange(plotlist = curr_plot, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1.5, ncol = 2, nrow =2, legend = "right", common.legend = TRUE)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS19.png", sep = ""), width = 800, height =600)

