##################################################
#  Preliminary Growth Curve
##################################################
## Project: FFAR Project
## Script purpose: Growth curves of selected preliminary isolates
##################################################
#Created by SL on 6-26-19
#Edited by
##################################################

growth_curve <- read.csv("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Growth curves/R Code/Practice_preliminary_growth_curves_SL_062619.csv")
library(ggplot2)
ggplot(data=growth_curve, aes(x=Hour, y=Concentration))+geom_point(aes(col=Isolate))
if(!require("devtools")) install.packages("devtools")
devtools::install_github("briandconnelly/growthcurve", build_vignettes = TRUE)
library(dplyr)
library(growthcurve)
rep1 <- filter(growth_curve, Strain=="FSL R10-0084")
myfit <- fit_growth(rep1, Hour, Concentration)
plot(myfit)
plot1<- ggplot(data=growth_curve, aes(x=Hour, y=Concentration, color=Isolate))+geom_point(shape=1)+xlab("Hours Post Inoculation")+ylab("CFU/mL")+ggtitle("Growth Curve of FFAR Strains")+stat_growthcurve()

#log transform
plot1 + scale_y_continuous(trans='log10')                                                                                                        
                                                                                                                  


#making individual plots
plot2<- ggplot(data=growth_curve, aes(x=Hour, y =Concentration)) + geom_point() +
  facet_wrap(~Isolate)
plot2+ scale_y_continuous(trans='log10')     

