## Figure S4: Final parametric model of prokaryotic abundance

library(dplyr)
library(ggplot2)
library(ggpubr)
library(caret)
library(reshape2)
library(visreg)

df <-read.csv('./data/curr_bac_final.csv', header=TRUE)
df <- df %>% dplyr::filter(silicate > 0 & phosphate > 0 & nitrate > 0 & oxygen > 0)

df$log_phosphate <- log10(df$phosphate)
df$log_nitrate <- log10(df$nitrate)
df$log_silicate <- log10(df$silicate)

df1 <- df %>% 
  dplyr::mutate(si_star = silicate - nitrate,
                n_star = nitrate-16*phosphate) %>%
  dplyr::select(c("logabund","logchlo", "logdepth","temperature", "oxygen", "log_silicate", "log_nitrate", "log_phosphate", "n_star", "si_star", "aou"))%>%
  dplyr::filter(n_star < 50)


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
