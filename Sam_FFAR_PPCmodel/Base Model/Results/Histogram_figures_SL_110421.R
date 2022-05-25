#LOAD in RDATA file

library(readr)
library(MASS)
library(dplyr)
library(ggplot2)

#load data
load("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Initial Microbial Count/Analysis/CheckingInitialDistribution/Model with 1.11log initial count/Results/BASEMODEL_ppcmodel_1.11initial.RData")

day7<- data3 %>% filter(day == 7)
day10<- data3 %>% filter(day == 10)
day14<- data3 %>% filter(day == 14)



#attempt to plot out initial microbial count data
#histogram with frequencies
jpeg(file="day7histogram.jpeg")
hist(day7$count, 10, main="Day 7", 
     xlab=(expression("LOG"[10]~"CFU/mL")), xlim=c(-5,10), ylim = c(0,15000), freq = TRUE, breaks = 20)
dev.off()

hist(day10$count, 10, main="Day 10",
     xlab=(expression("LOG"[10]~"CFU/mL")), xlim=c(-5,10), ylim = c(0,15000), freq = TRUE, breaks = 20)

hist(day14$count, 10, main="Day 14",
     xlab=(expression("LOG"[10]~"CFU/mL")), xlim=c(-5,10), ylim = c(0,15000), freq = TRUE, breaks = 20)

#attempt to put all histograms in one place
day7count<- as.data.frame(day7$count)
day10count<- as.data.frame(day10$count)
day14count<- as.data.frame(day14$count)

colnames(day7count)[1] <- "count"
colnames(day10count)[1] <- "count"
colnames(day14count)[1] <- "count"
day7count$day<-"7"
day10count$day<-"10"
day14count$day<-"14"

daycounts<-rbind(day7count,day10count,day14count)

ggplot(daycounts, aes(count, fill = day)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..),position = 'identity')