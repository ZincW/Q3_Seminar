library(foreign)
library(car) #Package includes Levene's test 
library(tidyr)  # for wide to long format transformation of the data
library(ggplot2)
library(QuantPsyc) #include lm.beta()
library(gmodels)
library(MASS)
sink("Part2Output.txt")
SalesData <- read.csv("/Users/zina/Google Drive/Q3_2017/Q3 Seminar/saledata.csv")
colnames(SalesData) <- c("RetailSales"," RetailEstablishments"," Income"," FederalExpenditures"," GenderRatio")
SalesData = SalesData[,c(2,3,4,5)]
attach(SalesData)
################ distribution of dependent variable #########
library(ggplot2)
ggplot(SalesData,aes(x=RetailEstablishments))+geom_histogram()
hist(RetailEstablishments)
plot(density(RetailEstablishments))

par(mfrow=c(2,2))
plot(Income,RetailEstablishments)
plot(FederalExpenditures,RetailEstablishments)
plot(GenderRatio,RetailEstablishments)

################# Chi-square #########
CrossTable(Gender, HeadWeight, fisher = FALSE, chisq = TRUE, expected = TRUE, sresid = TRUE)
#################### multiple linear regression ####################

SaleModel<-lm(RetailEstablishments~Income+FederalExpenditures+GenderRatio,data=SalesData)
summary(SaleModel)
confint(SaleModel,level=0.95)
coef_lmbeta <- lm.beta(SaleModel)
coef_lmbeta
##################### Examine Assumption ####################
tapply(HeadWeight,HeadSize, shapiro.test) # test normaliy of each level
leveneTest(RetailEstablishments)
plot(SaleModel)
shapiro.test(HeadWeight)
#testing independence of error
durbinWatsonTest(SaleModel)
#test of multicollinearity
vif(SaleModel)
1/vif(SaleModel) # Tolerance
### histogram of studentized residual to check if normal distributed
par(mfrow=c(1,2))
hist(SaleModel$residuals)
hist(rstudent(SaleModel))
plot(SaleModel$residuals, SaleModel$fitted)
plot(SaleModel)
sink()