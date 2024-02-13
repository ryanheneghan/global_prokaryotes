
library(dplyr)
library(hacksaw)
library(sp)
library(stringr)
library(mgcv)
library(visreg)
library(lme4)
library(MuMIn)


### BUILD STATISTICAL MODEL FOR SGRs
all_data <- read.csv("./data/all_sgr.csv")
all_data <- all_data %>% dplyr::mutate(across(-c(source), as.numeric))
all_data$logSGR <- log10(all_data$SGR)
all_data$logCHL <- log10(all_data$CHL)

gam_sst <- gam(logSGR ~ s(SST) + source, data = all_data)
gam_chl <- gam(logSGR ~ s(SST) + s(logCHL) + source, data = all_data)
par(mfrow = c(2,1))
visreg(gam_sst, xvar = "SST", residuals = TRUE)
visreg(gam_chl, xvar = "logCHL", residuals = TRUE)
save(gam_sst, file = "./models/gam_sst_sgr.RData")
save(gam_chl, file = "./models/gam_chl_sgr.RData")

lmer_sst <- lmer(logSGR ~ SST + (1|source), data = all_data) 
lmer_chl <- lmer(logSGR ~ SST + logCHL + (1|source), data = all_data) 
par(mfrow = c(2,1))
visreg(lmer_sst, xvar = "SST", ylab = "Specific Growth Rate (d^-1)", xlab = "Temperature (C)")
visreg(lmer_chl, xvar = "logCHL", ylab = "Specific Growth Rate (d^-1)", xlab = "log10(Chlo, mg m3)")
save(lmer_sst, file = "./models/lmer_sst_sgr.RData")
save(lmer_chl, file = "./models/lmer_chl_sgr.RData")


#### Surface waters, random effect for data source
surf_data <- all_data %>% filter(!is.na(CHL))
lmer_surf <- lmer(logSGR ~ SST + logCHL + (1|source), data = surf_data) 
summary(lmer_surf)
par(mfrow = c(1,1))
visreg(lmer_surf, "SST")
visreg(lmer_surf, "logCHL", ylab = "Specific Growth Rate (d^-1)", xlab = "log10(Chlo, mg m3)")


#### All waters, temperature only, random effect for data source
deep_data <- all_data %>% filter(is.na(CHL))
lmer_deep <- lmer(logSGR ~ SST + (1|source), data = deep_data) 
summary(lmer_deep)
par(mfrow = c(1,1))
visreg(lmer_deep, "SST")

par(mfrow = c(2,1))
visreg(lmer_surf, "SST")
visreg(lmer_deep, "SST")


