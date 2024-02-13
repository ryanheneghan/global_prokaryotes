
library(dplyr)
library(tidyr)
library(caret)
library(randomForest)
library(mgcv)
library(visreg)
library(car)
library(lme4)
library(segmented)
library("MuMIn")


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
                    Leg_5 = ifelse(Leg == "5", "Yes", "No"),
                    log10_depth = ifelse(log10_depth < 1, 1, log10_depth)) %>% 
        drop_na(Leg, Temp_C, log10_depth, log10_BT)#%>% filter(Leg != "5")

#Split the data set, sample, into a training data set (80%) and test data set (20%).
#The training data is used to train the model, and then the model is run on the test data to check how well it is classifying the data.

inTrain <- createDataPartition(y=df$fgC_cell,
                               times = 1,
                               list = FALSE,
                               p = .8)

training <- df[inTrain,]
testing <- df[-inTrain,]

#Fit random forest for variable selection

# make predictions
x_test <- dplyr::select(testing, -fgC_cell)
head(x_test) #x_test now has 3 variables, doesn't include Class
y_test <- dplyr::select(testing, fgC_cell) #Choose the column for Class


### GAM TO EXPLORE RELATIONSHIP
gam1 <-gam(fgC_cell ~ s(Temp_C) + Station, data = training)
summary(gam1)
par(mfrow = c(1,2))
visreg(gam1)
save(gam1, file = "./models/gam_model1_carbon.RData")

### GLMM WITH TEMPERATURE AND RANDOM EFFECT FOR STATION
glmm1 <- lmer(fgC_cell ~ poly(Temp_C,3) +(1|Station), data = training)
summary(glmm1)
r.squaredGLMM(glmm1)
visreg(glmm1)
save(glmm1, file = "./models/glmm_model1_carbon.RData")


### GLMM, BUT WITH STATION 5 REMOVED
df_sub <- df %>% dplyr::filter(Leg != "5")
glmm2 <- lmer(fgC_cell ~ poly(Temp_C,3)+(1|Station), data = df_sub)
summary(glmm2)
r.squaredGLMM(glmm2)
visreg(glmm2)

### GLM, NO RANDOM EFFECT FOR STATION
glmm2 <- lm(fgC_cell ~ poly(Temp_C,3)+poly(log10_depth,3), data = training)
check_collinearity(glmm2)
summary(glmm2)
#r.squaredGLMM(glmm2)
visreg(glmm2)
vif(glmm2)
