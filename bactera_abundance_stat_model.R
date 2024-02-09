rm(list = ls())

#Example random forest variable selection - J Holloway-Brown August 2021
library(dplyr)
library(caret)
library(randomForest)
library(mgcv)
library(car)

setwd("~/Desktop/Papers/Bacteria_Census")

df <-read.csv('./curr_bac_final.csv', header=TRUE)
df <- df %>% dplyr::filter(silicate > 0 & phosphate > 0 & nitrate > 0 & oxygen > 0)

#df[df$phosphate < 10^-3, "phosphate"] <- 10^-3
df$log_phosphate <- log10(df$phosphate)

#df[df$nitrate < 10^-4, "nitrate"] <- 10^-4
df$log_nitrate <- log10(df$nitrate)

#df[df$silicate < 10^-4, "silicate"] <- 10^-4
df$log_silicate <- log10(df$silicate)

df1 <- df %>% #dplyr::filter(10^logdepth > z_eu) %>% 
dplyr::mutate(si_star = silicate - nitrate,
              n_star = nitrate-16*phosphate) %>%
  dplyr::select(c("logabund","logchlo", "logdepth","temperature", "oxygen", "log_silicate", "log_nitrate", "log_phosphate", "n_star", "si_star", "aou"))%>%
  dplyr::filter(n_star < 50)



# %>% dplyr::mutate(n_star = nitrate - 16*phosphate)
#df1 <- df %>% dplyr::filter(10^logdepth >= z_eu) %>% dplyr::select(c("logabund","logdepth", "temperature", "log_s_phosphate", "log_s_nitrate", "oxygen", "silicate"))
#df1[df1$phosphate < 1e-2,"phosphate"] <- 1e-2
#df1[df1$nitrate < 1e-2, "nitrate"] <- 1e-2
#df1[df1$silicate < 1e-2, "silicate"] <- 1e-2

#df1 <- df1 %>% mutate(phosphate = log10(phosphate),
#                      nitrate = log10(nitrate),
#                      silicate = log10(silicate))

#Split the data set, sample, into a training data set (80%) and test data set (20%).
#The training data is used to train the model, and then the model is run on the test data to check how well it is classifying the data.

inTrain <- createDataPartition(y=df1$logabund,
                               times = 1,
                               list = FALSE,
                               p = .8)

training <- df1 [inTrain,]
testing <- df1 [-inTrain,]


#write.csv(training, file= "training_bacteria.csv")
#write.csv(testing, file="testing_bacteria.csv")

#Fit random forest for variable selection

# make predictions
x_test <- dplyr::select(testing, -logabund)
head(x_test) #x_test now has 3 variables, doesn't include Class
y_test <- dplyr::select(testing, logabund) #Choose the column for Class

set.seed(1) #you can choose any seed but the structure of the trees will be slightly different each time if you don't set one

#1 specify the random forest model for all variables
rf <- randomForest(logabund ~., data=training, importance = TRUE)

write.csv(rf$importance, "node_purity.csv", row.names = FALSE)

# 2 - Run predictions and return results as percentages
rf.pred <- predict(rf, x_test, type = "response")

#plot variable importance - most important variable to least important variable from top right to bottom left of plot. Can also view the results as a table
par(mfrow = c(1,1))
varImpPlot(rf, sort=TRUE)

rf.pred <- predict(rf, x_test, type = "response")

#Produce accuracy metrics for the random forest- we want models with low error because if a model has very poor accuracy, the variable importance may not be reliable

error<- (rf.pred - y_test)
error<-as.numeric(error[,1])

#Function to calculate mean square error
rmse <- function(error)
{
  sqrt(mean(error)^2)
}

rmse(error)

#Calculate mean absolute error
mae <- function(error)
{
  mean(abs(error))
}

mae(error)


save(rf,file = "rf_model2.RData")

#I usually produce multiple random forests with different seeds- at least 5 or enough to show there is a relatively consistent pattern to the variable importance

#Next step is to choose the variables with the highest importance and proceed to parametric modelling

gam1 <- gam(logabund ~ s(logdepth)+s(log_nitrate)+
              s(aou) + s(temperature) +
              s(logchlo), data = training)
concurvity(gam1)

summary(gam1)


rf.pred <- predict(gam1, x_test, type = "response")

error<- (rf.pred - y_test)
error<-as.numeric(error[,1])

rmse(error)
mae(error)



par(mfrow = c(3,2))
visreg(gam1)

#par(mfrow = c(1,1))
#visreg2d(gam1,"logdepth","aou",plot.type="persp")
#visreg2d(gam1,"logdepth","log_nitrate",plot.type="persp")
#visreg2d(gam1,"logdepth","temperature",plot.type="persp")
#visreg2d(gam1,"logdepth","logchlo",plot.type="persp")

save(gam1, file = "gam_model_final.RData")




## FIT PARAMETRIC MODEL
gm1 <- lm(logabund ~ poly(logdepth,5) + log_nitrate + poly(aou, 2) + logchlo + temperature, data = training)
summary(gm1)
vif(gm1)

save(gm1, file = "lm_model_final.RData")

par(mfrow = c(3,2))
visreg(gm1)
par(mfrow = c(1,1))
varImpPlot(rf)

rf.pred <- predict(gm1, x_test, type = "response")

error<- (rf.pred - y_test)
error<-as.numeric(error[,1])

rmse(error)
mae(error)


## FIT SURFACE AND DEEP WATER MODELS SEPERATELY
training <- training %>% dplyr::select(logabund, log_nitrate, aou, logchlo, temperature, logdepth)
training <- training %>% dplyr::mutate(logdepth2 = logdepth^2, logdepth3 = logdepth^3, logdepth4 = logdepth^4, logdepth5 = logdepth^5,
                                       aou2 = aou^2)
training_surface <- training %>% dplyr::filter(logdepth <= log10(100)) 
training_deep <- training %>% dplyr::filter(logdepth > log10(100)) 

lm_surface <- lm(logabund ~ poly(logdepth, 1) + log_nitrate+ poly(aou, 2) +  temperature + logchlo, data = training_surface)
summary(lm_surface)
par(mfrow = c(2,3))
visreg(lm_surface, ylim = c(4.5,7))

lm_deep <- lm(logabund ~logdepth + logdepth2 + logdepth3 + logdepth4+ logdepth5+ log_nitrate +  aou + aou2+ logchlo + temperature, data = training_deep)
summary(lm_deep)
par(mfrow = c(3,2))
visreg(lm_deep)

lm_all <- lm(logabund ~ logdepth + logdepth2 + logdepth3 + logdepth4+ logdepth5+ log_nitrate +  aou + aou2+ logchlo + temperature, data = training)
summary(lm_all)
par(mfrow = c(3,2))
visreg(lm_all)

par(mfrow = c(1,1))

cor(training_epi)
cor(training_deep)
