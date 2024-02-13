### Figure S1: RAW RELATIONSHIP BETWEEN BACT ABUNDANCE AND FULL PREDICTOR VARIABLES
library(viridis)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)
library(dplyr)
library(caret)
library(reshape2)

df <-read.csv('./data/curr_bac_final.csv', header=TRUE)

df <- df %>% dplyr::filter(silicate > 0 & phosphate > 0 & nitrate > 0 & oxygen > 0)

df$log_phosphate <- log10(df$phosphate)
df$log_nitrate <- log10(df$nitrate)
df$log_silicate <- log10(df$silicate)

df1 <- df %>% 
  dplyr::mutate(si_star = silicate - nitrate,
                n_star = nitrate-16*phosphate) %>%
  dplyr::select(c("logabund", "logdepth","logchlo","temperature", "oxygen", "log_silicate", "log_nitrate", "log_phosphate", "n_star", "si_star", "aou"))%>%
  dplyr::filter(n_star < 50)

x_names <- c("logdepth","logchlo", "temperature", "oxygen",  "log_nitrate", "log_phosphate", "log_silicate", "aou", "n_star", "si_star")
x_save_name <- c(expression(bold(paste("log"[10], "(Sample Depth, m)", sep = ""))),
                 expression(bold(paste("log"[10], "(Chlorophyll mg m"^-3, ")", sep = ""))),
                 expression(bold("Temperature, °C" )),
                 expression(bold(paste("Oxygen (", mu, "mol kg"^-1, ")", sep = ""))),
                 expression(bold(paste("log"[10], "(Nitrate ", mu, "mol kg"^-1, ")", sep = ""))),
                 expression(bold(paste("log"[10], "(Phosphate ", mu, "mol kg"^-1, ")", sep = ""))),
                 expression(bold(paste("log"[10], "(Silicate ", mu, "mol kg"^-1, ")", sep = ""))),
                 expression(bold(paste("AOU (", mu, "mol kg"^-1, ")", sep = ""))),
                 expression(bold("Si*")),
                 expression(bold("N*")))

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
  curr_dat <- df1[,c("logabund",x_names[i])]
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
    ylab(expression(bold(paste("log"[10], "(Bacteria, ml"^-1, ")", sep = ""))))+ 
    xlab(x_save_name[i]) + guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5)) +
    theme_opts
}



tagss <-c("a", "b", "c", "d", "e", "f","g","h","i","j")
figure <- ggpubr::ggarrange(plotlist = curr_plot, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1.5, ncol = 3, nrow = 4, legend = "right", common.legend = TRUE)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS1.png", sep = ""), width = 1250, height = 1470)

