#### Figure S15: REGIONAL ABUNDANCE MODEL ASSESSMENT

library(dplyr)
library(caret)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(visreg)
library(xlsx)

df <- read.csv('./data/raw_observations/prokaryote_abundance.csv') # Import prokaryote abundance data
df <- df %>% dplyr::filter(silicate > 0 & phosphate > 0 & nitrate > 0 & oxygen > 0)

df$log_phosphate <- log10(df$phosphate)
df$log_nitrate <- log10(df$nitrate)
df$log_silicate <- log10(df$silicate)

df1 <- df %>% 
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

tagss <-c("a", "b", "c", "d")

plot_list <- list()
wb <- createWorkbook(type = "xls")

for(i in 1:length(regions)){
  reg_dat1 <- all1 %>% dplyr::filter(abs(lat) > min_lat[i] & abs(lat) <= max_lat[i])
  x_test1 <- reg_dat1 %>% dplyr::select(logdepth,logchlo, temperature, log_nitrate,aou)
  reg_pred1 <- predict(gm1, x_test1, type = "response")
  error1<- (reg_pred1 - reg_dat1$logabund)

  gg_dat <- data.frame("x_dat" = reg_dat1$logabund, "y_dat" = reg_pred1)
  curr_r2 <- cor(gg_dat$x_dat, gg_dat$y_dat)^2
  ann_label <- expression("R"^2 == curr_r2)
  
  plot_list[[i]] <- ggplot(gg_dat) + geom_point(aes(x=x_dat,y=y_dat), shape = 1) + geom_abline(slope = 1, intercept = 0, colour = "red", size = 1) +
    theme_bw() + theme_opts + xlab(expression(bold(paste("Observed log"[10], "(abundance, mL"^-1, ")", sep = "")))) +
    ylab(expression(bold(paste("Predicted log"[10], "(abundance, mL"^-1, ")", sep = "")))) +
    labs(title = paste(regions[i], ", n = ", length(reg_pred1), sep = "")) + xlim(c(min(df1$logabund),max(df1$logabund))) + ylim(c(4.3,6.2))+
    annotate("label", x =6.5, y = 4.4, parse = TRUE, label =paste("R ^ 2 == ", round(curr_r2*100,1),sep = ""), size = 6, fontface = "bold")
  
  ## Save source data for Figure S15
  sheet_data <- gg_dat %>% dplyr::rename(Observed_log10_prokaryote_abundance_mL = x_dat, Predicted_log10_prokaryote_abundance_mL = y_dat)
  header_name <- paste("Figure S15", tagss[i], sep = "")
    
  sh1 <- createSheet(wb,sheetName = header_name)
  addDataFrame(data.frame(header_name=c(header_name)),sheet = sh1,row.names = FALSE,startRow = 1,col.names = FALSE)
  addDataFrame(sheet_data,sheet = sh1,row.names = FALSE,startRow = 2)
}

saveWorkbook(wb,"./source_data/FigureS15.xls")

figure <- ggpubr::ggarrange(plotlist = plot_list, labels = tagss, font.label = list(size = 22), hjust = -0.5, 
                            vjust = 1 , ncol = 2, nrow = 2)
ggpubr::ggexport(figure, filename = paste("./figures/FigureS15.png", sep = ""), width = 800, height = 700)

