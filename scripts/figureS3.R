### FIGURE S3: Random forest variable importance plot

library(dplyr)
library(ggplot2)
library(ggpubr)
library(caret)
library(randomForest)
library(mgcv)
library(reshape2)
library(visreg)


theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "black", linewidth=0.7),
                         axis.text = element_text(size = 14, face = "bold", colour = "black"),
                         axis.title = element_text(size = 16, face = "bold", colour = "black"),
                         plot.title = element_text(size = 18, face = "bold.italic", colour = "black"),
                         plot.margin = unit(c(0.6,0.3,0.3,0.3), "cm")))

## VARIABLE IMPORTANCE PLOT
imp_data <- read.csv("./data/node_purity.csv")
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


ggpubr::ggexport(figure_s2, filename = paste("./figures/FigureS3.png", sep = ""), width = 600, height = 400)

## Save source data
df_original <-read.csv('./data/node_purity.csv', header=TRUE)
write.csv(df_original, "./source_data/FigureS3/FigureS3.csv", row.names = FALSE)
