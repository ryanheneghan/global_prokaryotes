### Figure S9: GLM FOR PROK CELL CARBON (FIXED EFFECT ONLY)

library(dplyr)
library(tidyr)
library(viridis)
library(ggpointdensity)
library(ggplot2)
library(ggpubr)
library(lme4)

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

inTrain <- createDataPartition(y=df$fgC_cell,
                               times = 1,
                               list = FALSE,
                               p = .8)

training <- df[inTrain,]

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

gm1 <- lmer(fgC_cell ~ poly(Temp_C,3)+(1|Station), data = training)
summary(gm1)

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

