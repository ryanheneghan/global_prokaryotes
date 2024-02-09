rm(list = ls())

library(dplyr)
library(tidyr)
library(viridis)
library(ggpointdensity)
library(ggplot2)
library(ggpubr)
library(lme4)

setwd("~/Desktop/Papers/Bacteria_Census")

df <-read.csv('malaspina_individual_bacteria.csv', header=TRUE)
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

### RAW RELATIONSHIP BETWEEN BACT ABUNDANCE AND FULL PREDICTOR VARIABLES


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
dot_dat <- data.frame("Temp_C" = tt$res$Temp_C, "Value" = tt$res$visregRes)

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
dot_dat <- data.frame("Station" = tt$res$Station, "Value" = tt$res$visregRes)

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


#### FIGURE OF FINAL PARAMETRIC MODEL
gm1 <- lmer(fgC_cell ~ poly(Temp_C,3)+(1|Station), data = training)
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
dot_dat <- data.frame("Temp_C" = tt$res$Temp_C, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$Temp_C,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[1]] <- ggplot() + geom_point(data=dot_dat, aes(x=Temp_C, y = Value, colour = dens), size = 1)+ 
  geom_line(data=line_dat, aes(x=Temp_C, y = Value), col = "red", size = 2) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = Temp_C, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold("Temperature, °C" )))+
  ylab(expression(bold("fg C cell"^{-1}))) + 
  theme_opts


figure <- ggpubr::ggarrange(plotlist = plot_list,  ncol = 1, nrow = 1, legend = "right", common.legend = TRUE)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS9.png", sep = ""), width = 500, height = 400)



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

