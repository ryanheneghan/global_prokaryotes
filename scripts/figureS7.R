### Figure S7: RAW RELATIONSHIP BETWEEN PROK CELL CARBON AND FULL PREDICTOR VARIABLES

library(dplyr)
library(tidyr)
library(viridis)
library(ggpointdensity)
library(ggplot2)
library(ggpubr)

df <-read.csv('./data/malaspina_individual_bacteria.csv', header=TRUE)
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


x_names <- c("log10_depth","Temp_C", "Station")
x_save_name <- c(expression(bold(paste("log"[10], "(Sample Depth, m)", sep = ""))),
                 expression(bold("Temperature, °C" )),
                 expression(bold("Station")))

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

curr_plot <- list()


for(i in 1:length(x_names)){
  curr_dat <- df[,c("fgC_cell",x_names[i])]
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
  
  if(x_names[i] == "Station"){
    curr_plot[[i]] <- curr_plot[[i]] + scale_x_discrete(x_save_name[i], breaks = c("001", "070", "140"), labels = c("1", "70", "140"))
  }
}



tagss <-c("a", "b", "c")
figure <- ggpubr::ggarrange(plotlist = curr_plot, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1.5, ncol = 2, nrow =2, legend = "right", common.legend = TRUE)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS7.png", sep = ""), width = 800, height =600)

