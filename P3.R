library(foreign)
library(car)
library(ggplot2)
library(nlme)
library(reshape)
library(graphics)
library(readr)
library(readr)
sink("Part3Output.txt")
set2 <- read_csv("Downloads/set2.csv")
colnames(set2) <- c("Subject", "Session",'Score')
attach(set2)
############## Distribution #########
par(mfrow=c(2,2))
hist(Score,main="Histogram of Score", xlab="Score",data=set2)
plot(density(Score),main="Density of Score",data=set2)
plot(Score,main="qqPlot of Score",data=set2)
boxplot(Score,main="Boxplot of Score",data=set2)
############## Relationship ###########
plot(Session,Score)
abline(lm(Score ~ Session))
############# Random Intercept #########
randomInterceptOnly <- lme(Score ~ 1, data = set2, random = ~1|Subject, method = "ML")

############ Random Intercept and Fixed Factors ############
RandomInterceptFixedFactors <- lme(Score ~ Session, data = set2, random = ~1|Subject, method = "ML")

############ Random slopes model ###############
RandomSlopeIntercept <-lme(Score ~ Session, data = set2, random = ~Session|Subject, method = "ML",control=lmeControl(opt = "optim") )

summary(randomInterceptOnly)
summary(RandomInterceptFixedFactors)
summary(RandomSlopeIntercept)
anova(randomInterceptOnly,RandomInterceptFixedFactors,RandomSlopeIntercept)
intervals(randomInterceptOnly, 0.95) 
summary(randomInterceptOnly)

CheckVariance <-lme(Score ~ Subject, data = set2, random = ~Session|Subject, method = "ML",control=lmeControl(opt = "optim") )
intervals(CheckVariance, 0.95) 
summary(CheckVariance)
sink()