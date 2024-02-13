#### FIGURE S5: Sensitivity analysis of abundance model

library(ggplot2)
library(ggpubr)
library(dplyr)
library(caret)
library(reshape2)
library(visreg)


df <-read.csv('./data/curr_bac_final.csv', header=TRUE)
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
