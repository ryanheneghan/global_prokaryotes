rm(list = ls())

library(dplyr)
library(caret)
library(randomForest)
library(mgcv)
library(reshape2)
library(visreg)

setwd("~/Desktop/Papers/Bacteria_Census")

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", linewidth=0.7),
                         axis.text = element_text(size = 14, face = "bold", colour = "black"),
                         axis.title = element_text(size = 16, face = "bold", colour = "black"),
                         plot.title = element_text(size = 18, face = "bold.italic", colour = "black"),
                         plot.margin = unit(c(0.6,0.3,0.3,0.3), "cm")))

## VARIABLE IMPORTANCE PLOT
imp_data <- read.csv("./node_purity.csv")
imp_data <- imp_data[order(imp_data$IncNodePurity),]
imp_data$var_num <- 1:10

var_name <- c(expression(bold("Si*")), 
              expression(bold("Oxygen")), 
              expression(bold(paste("log"[10], "(Chlorophyll)", sep = ""))),
                 expression(bold("N*")), 
              expression(bold("Temperature")), 
              expression(bold(paste("log"[10], "(Silicate)", sep = ""))),
                expression(bold("AOU")), 
              expression(bold(paste("log"[10], "(Phosphate)", sep = ""))),
                   expression(bold(paste("log"[10], "(Nitrate)", sep = ""))),
                   expression(bold(paste("log"[10], "(Depth)", sep = ""))))

figure_s2 <- ggplot(imp_data) + geom_point(aes(x=IncNodePurity, y = var_num), size = 2) + theme_bw() +
  scale_y_continuous(name = "", breaks = 1:10, labels = var_name) + xlab("Increase in Node Purity") + theme_opts


ggpubr::ggexport(figure_s2, filename = paste("./figures/FigureS2.png", sep = ""), width = 600, height = 400)


#### FIGURE S5
df <-read.csv('./curr_bac_final.csv', header=TRUE)
df <- df %>% dplyr::filter(silicate > 0 & phosphate > 0 & nitrate > 0 & oxygen > 0)

df$log_phosphate <- log10(df$phosphate)
df$log_nitrate <- log10(df$nitrate)
df$log_silicate <- log10(df$silicate)

df1 <- df %>% #dplyr::filter(10^logdepth > z_eu) %>% 
  dplyr::mutate(si_star = silicate - nitrate,
                n_star = nitrate-16*phosphate) %>%
  dplyr::select(c("logabund", "logdepth","logchlo", "temperature", "oxygen", "log_silicate", "log_nitrate", "log_phosphate", "n_star", "si_star", "aou"))%>%
  dplyr::filter(n_star < 50)

inTrain <- createDataPartition(y=df1$logabund,
                               times = 1,
                               list = FALSE,
                               p = .8)

df1 <- df1 [inTrain,]

####### MODEL PARAMETER AND R2 SENSITIVITY TO SAMPLE SIZE
num_vals <- c(1000,5000,10000,20000,30000)
num_samp <- 1000
df1$logdepth2 <- (df1$logdepth)^2
df1$logdepth3 <- (df1$logdepth)^3
df1$logdepth4 <- (df1$logdepth)^4
df1$logdepth5 <- (df1$logdepth)^5
df1$aou2 <- (df1$aou)^2
save_out <- array(NA, dim = c(num_samp, length(num_vals), 12))
dimnames(save_out)[[2]] <- c("1,000", "5,000", "10,000", "20,000", "30,000")
dimnames(save_out)[[3]] <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k","R^2")

for(i in 1:length(num_vals)){
  curr_num <- num_vals[i]
  print(curr_num)
  
  pb = txtProgressBar(min = 0, max =num_samp, initial = 0, style = 3)
  
  for(j in 1:num_samp){
    setTxtProgressBar(pb,j)
    sample2 <- df1[sample(1:nrow(df1), curr_num,
                          replace=TRUE),]
    gm1 <- lm(logabund ~ logdepth+ logdepth2 + logdepth3 + logdepth4 + logdepth5 +
                          log_nitrate + aou + aou2 + temperature + logchlo, data = sample2)
    kk <- summary(gm1)
    
    save_out[j,i,1] <- gm1$coefficients[1]
    save_out[j,i,2] <- gm1$coefficients[2]
    save_out[j,i,3] <- gm1$coefficients[3]
    save_out[j,i,4] <- gm1$coefficients[4]
    save_out[j,i,5] <- gm1$coefficients[5]
    save_out[j,i,6] <- gm1$coefficients[6]
    save_out[j,i,7] <- gm1$coefficients[7]
    save_out[j,i,8] <- gm1$coefficients[8]
    save_out[j,i,9] <- gm1$coefficients[9]
    save_out[j,i,10] <- gm1$coefficients[10]
    save_out[j,i,11] <- gm1$coefficients[11]
    save_out[j,i,12] <- kk$adj.r.squared
  }
}





theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", linewidth=0.7),
                         axis.text = element_text(size = 14, face = "bold", colour = "black"),
                         axis.title = element_text(size = 16, face = "bold", colour = "black"),
                         plot.title = element_text(size = 18, face = "bold.italic", colour = "black"),
                         plot.margin = unit(c(0.6,0.3,0.3,0.3), "cm")))

curr_plot <- list()

for(i in 1:length(dimnames(save_out)[[3]])){
  gg_curr_dat <- save_out[,,i]
  gg_curr_dat <- melt(gg_curr_dat)
  colnames(gg_curr_dat) <- c("Var1", "Num_Obs", "Value")
  curr_name <- dimnames(save_out)[[3]][i]
  
  curr_plot[[i]] <- ggplot(gg_curr_dat, aes(x=Num_Obs, y=Value)) + 
    geom_boxplot() + theme_bw() +
    xlab("") + ylab("") + labs(title = curr_name)+
    theme_opts
  
  if(curr_name == "R^2"){
    curr_plot[[i]] <- ggplot(gg_curr_dat, aes(x=Num_Obs, y=Value)) + 
      geom_boxplot() + theme_bw() +
      xlab("") + ylab("") + labs(title = expression(bolditalic("R"^2)))+
      theme_opts
  }
  
}


tagss <-c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k","l")
figure <- ggpubr::ggarrange(plotlist = curr_plot, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1.5, ncol = 3, nrow = 4)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS5.png", sep = ""), width = 1200, height = 1480)


### RAW RELATIONSHIP BETWEEN BACT ABUNDANCE AND FULL PREDICTOR VARIABLES
library(viridis)
library(ggpointdensity)

df <-read.csv('./curr_bac_final.csv', header=TRUE)
df <- df %>% dplyr::filter(silicate > 0 & phosphate > 0 & nitrate > 0 & oxygen > 0)

df$log_phosphate <- log10(df$phosphate)
df$log_nitrate <- log10(df$nitrate)
df$log_silicate <- log10(df$silicate)

df1 <- df %>% #dplyr::filter(10^logdepth > z_eu) %>% 
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


#### FIGURE OF GAM RELATIONSHIPS
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

gm1 <- gam(logabund ~  s(logdepth)+s(log_nitrate)+
             s(aou) + s(temperature) +
             s(logchlo), data = df1)
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

tt = visreg(gm1, xvar = "logdepth", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("logdepth" = tt$fit$logdepth, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("logdepth" = tt$res$logdepth, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$logdepth,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[1]] <- ggplot() + geom_point(data=dot_dat, aes(x=logdepth, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=logdepth, y = Value), col = "red", size = 1.5) + 
  geom_ribbon(data = line_dat, aes(x = logdepth, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  theme_bw()+
  xlab(expression(bold(paste("log"[10], "(Depth, m)", sep = "")))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  theme_opts

tt = visreg(gm1, xvar = "log_nitrate", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("log_nitrate" = tt$fit$log_nitrate, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("log_nitrate" = tt$res$log_nitrate, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$log_nitrate,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[2]] <- ggplot() + geom_point(data=dot_dat, aes(x=log_nitrate, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=log_nitrate, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = log_nitrate, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold(paste("log"[10],"(Nitrate ", mu, "mol kg"^-1, ")", sep = "")))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+
  theme_opts

tt = visreg(gm1, xvar = "aou", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("aou" = tt$fit$aou, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("aou" = tt$res$aou, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$aou,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[3]] <- ggplot() + geom_point(data=dot_dat, aes(x=aou, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=aou, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = aou, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold(paste("AOU (", mu, "mol kg"^-1, ")", sep = "")))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+
  theme_opts

tt = visreg(gm1, xvar = "temperature", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("temperature" = tt$fit$temperature, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("temperature" = tt$res$temperature, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$temperature ,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[4]] <- ggplot() + geom_point(data=dot_dat, aes(x=temperature, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=temperature, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = temperature, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold("Temperature, °C" ))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+
  theme_opts

tt = visreg(gm1, xvar = "logchlo", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("logchlo" = tt$fit$logchlo, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("logchlo" = tt$res$logchlo, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$logchlo ,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[5]] <- ggplot() + geom_point(data=dot_dat, aes(x=logchlo, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=logchlo, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = logchlo, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold(paste("log"[10], "(Chlorophyll, mg m"^-3, ")", sep = "")))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+
  theme_opts

tagss <-c("a", "b", "c", "d", "e")
figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1 , ncol = 3, nrow = 2, legend = "right", common.legend = TRUE)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS3.png", sep = ""), width = 1250, height = 700)


#### FIGURE OF FINAL PARAMETRIC MODEL
gm1 <- lm(logabund ~ poly(logdepth,5) + log_nitrate + poly(aou, 2) + temperature + logchlo, data = df1)
summary(gm1)

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", size=0.7),
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

tt = visreg(gm1, xvar = "logdepth", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("logdepth" = tt$fit$logdepth, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("logdepth" = tt$res$logdepth, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$logdepth,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[1]] <- ggplot() + geom_point(data=dot_dat, aes(x=logdepth, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=logdepth, y = Value), col = "red", size = 1.5) + 
  geom_ribbon(data = line_dat, aes(x = logdepth, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  theme_bw()+
  xlab(expression(bold(paste("log"[10], "(Depth, m)", sep = "")))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  theme_opts

tt = visreg(gm1, xvar = "log_nitrate", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("log_nitrate" = tt$fit$log_nitrate, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("log_nitrate" = tt$res$log_nitrate, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$log_nitrate,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[2]] <- ggplot() + geom_point(data=dot_dat, aes(x=log_nitrate, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=log_nitrate, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = log_nitrate, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold(paste("log"[10],"(Nitrate ", mu, "mol kg"^-1, ")", sep = "")))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+
  theme_opts

tt = visreg(gm1, xvar = "aou", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("aou" = tt$fit$aou, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("aou" = tt$res$aou, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$aou,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[3]] <- ggplot() + geom_point(data=dot_dat, aes(x=aou, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=aou, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = aou, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold(paste("AOU (", mu, "mol kg"^-1, ")", sep = "")))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+
  theme_opts

tt = visreg(gm1, xvar = "temperature", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("temperature" = tt$fit$temperature, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("temperature" = tt$res$temperature, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$temperature ,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[4]] <- ggplot() + geom_point(data=dot_dat, aes(x=temperature, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=temperature, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = temperature, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold("Temperature, °C" ))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+
  theme_opts

tt = visreg(gm1, xvar = "logchlo", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = TRUE, rug = FALSE)
line_dat <- data.frame("logchlo" = tt$fit$logchlo, "Value" = tt$fit$visregFit, "Lower" = tt$fit$visregLwr, "Upper" = tt$fit$visregUpr)
dot_dat <- data.frame("logchlo" = tt$res$logchlo, "Value" = tt$res$visregRes)

x <- densCols(dot_dat$logchlo ,dot_dat$Value, colramp=colorRampPalette(c("black", "white"))) ## Use densCols() output to get density at each point
dot_dat$dens <- col2rgb(x)[1,] + 1L 
dot_dat <- dot_dat[order(dot_dat$dens),] ## Reorder rows so that densest points are plotted on top

plot_list[[5]] <- ggplot() + geom_point(data=dot_dat, aes(x=logchlo, y = Value, colour = dens), size = 0.7)+ 
  geom_line(data=line_dat, aes(x=logchlo, y = Value), col = "red", size = 1.5) + theme_bw()+
  geom_ribbon(data = line_dat, aes(x = logchlo, ymin = Lower, ymax = Upper), fill = "red", alpha = .2) +
  scale_colour_viridis(name = "Num.\nSamples",  breaks = c(1,83,167,250), labels = c("1", "80",  "160",  "≥240"))+
  xlab(expression(bold(paste("log"[10], "(Chlorophyll, mg m"^-3, ")", sep = "")))) + 
  ylab(expression(bold(paste("log"[10], "(Abundance, ml"^-1, ")", sep = "")))) +
  ylim(c(4.5,6.5))+
  theme_opts

tagss <-c("a", "b", "c", "d", "e")
figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1 , ncol = 3, nrow = 2, legend = "right", common.legend = TRUE)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS4.png", sep = ""), width = 1250, height = 700)


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


tagss <-c("a", "b", "c", "d")
figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1 , ncol = 2, nrow = 2)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS6.png", sep = ""), width = 800, height = 700)

