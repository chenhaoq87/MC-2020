#Confidence intervals
#by SL on 9/16/19

#install.packages("minpack.lm")
#install.packages("nlstools")
library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
library(AICcmodavg)
library(minpack.lm)
library(nlstools)
#insert growth model data
#growthmodel2 <- read.csv("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Growth curves/R Code/Growth_Curve_Round4_RAWDATA_SL_090819.csv")
growthmodel2 <- read.csv("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Growth curves/R Code/Confidence intervals/Copy of Confidence_interval_RAWDATA_SL_091619.csv")
names(growthmodel2) <- c("t","Isolate","LOG10N","Count")

R101587 <- growthmodel2[growthmodel2$Isolate=="FSL R10-1587", ]

#Run Gompertz growth model
R101587modelg <- nlsLM(formula = gompertzm, 
                       data = R101587,
                       start = c(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                       lower=c(lag = 0, mumax=0, LOG10N0=0, LOG10Nmax=0))

#Get confidence interval
CI <- confint2(R101587modelg)

#You can view the results of this with
print(CI)

R100084 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0084", ]

#THIS DOESN'T WORK
#Run Gompertz growth model
#R100084modelg <- nls(gompertzm, R100084,
                     #list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))


R100084modelg <- nlsLM(formula = gompertzm, 
                       data = R100084,
                       start = c(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8), 
                       lower=c(lag = 0, mumax=0, LOG10N0=0, LOG10Nmax=0))



#Get confidence interval
CI2 <- confint2(R100084modelg)

#You can view the results of this with
print(CI2)

R103286 <- growthmodel2[growthmodel2$Isolate=="FSL R10-3286", ]

#Run Gompertz growth model
R103286modelg <- nlsLM(formula = gompertzm, 
                       data = R103286,
                       start = c(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                       lower=c(lag = 0, mumax=0, LOG10N0=0, LOG10Nmax=0))

#Get confidence interval
CI3 <- confint2(R103286modelg)

#You can view the results of this with
print(CI3)

R100941 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0941", ]

#Run Gompertz growth model
R100941modelg <- nlsLM(formula = gompertzm, 
                       data = R100941,
                       start = c(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                       lower=c(lag = 0, mumax=0, LOG10N0=0, LOG10Nmax=0))

#Get confidence interval
CI4 <- confint2(R100941modelg)

#You can view the results of this with
print(CI4)
