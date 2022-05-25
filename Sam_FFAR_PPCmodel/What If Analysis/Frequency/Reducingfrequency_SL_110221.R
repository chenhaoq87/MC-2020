load("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Initial Microbial Count/Analysis/CheckingInitialDistribution/What if analysis with 1.11log initial count/Frequency/BASEMODEL_ppcmodel_1.11initial.RData")

#make a duplicate of simulations
data4<-data3

n_units = n_sim*n_halfgal 

sampleKey <- data4[c(1,2)] %>%
  unique(.)


x= c(1,0)
prob100=c(1,0)
prob10 = c(0.1,0.9)
prob50 = c(0.5,0.5)
n= n_halfgal 
result100 = rep(x, round(n * prob100))
result10 = rep(x, round(n * prob10))
result50 = rep(x, round(n * prob50))

# sample n_units times at defined frequency of contamination
#(c,1,0) means it will assign one or zero to n units, n units is bulk tanks * half gallons which is equivalent 
#to simulations times half gallons, n unit= unique samples
#replace=true means you dont have a string of values

contam_BySample_freq100<-as.vector(replicate(n_sim, sample(result100)))
contam_BySample_freq50<-as.vector(replicate(n_sim, sample(result50)))
contam_BySample_freq10<-as.vector(replicate(n_sim, sample(result10)))

# bind your the string of values as a column to the sampleKey dataframe
cbind(sampleKey, contam_BySample_freq100) -> temp1
cbind(sampleKey, contam_BySample_freq50) -> temp2
cbind(sampleKey, contam_BySample_freq10) -> temp3

# merge your temporary dataframes with your full dataframe
merge(temp1,data4,by=c("BT_hg")) -> data5
merge(temp2,data5,by=c("BT_hg")) -> data5
merge(temp3,data5,by=c("BT_hg")) -> data5

# get updated count columns by multiplying value in contam_BySample_freqX by count
data5$count_freq100 <- data5$count*data5$contam_BySample_freq100 # should be the same as "count"
data5$count_freq50 <- data5$count*data5$contam_BySample_freq50 # should be the same as "count"
data5$count_freq10 <- data5$count*data5$contam_BySample_freq10 # should be the same as "count"


library(dplyr)
Day7Data <- data5[data5$day == "7", ]
Day10Data <- data5[data5$day == "10", ]
Day14Data <- data5[data5$day == "14", ]

#frequency reduced from 100% to 50%
sum(Day7Data$count_freq50 > 4.30103)
sum(Day10Data$count_freq50 > 4.30103)
sum(Day14Data$count_freq50 > 4.30103)

mean(Day7Data$count_freq50)
mean(Day10Data$count_freq50)
mean(Day14Data$count_freq50)

#frequency reduced from 100% to 10%
sum(Day7Data$count_freq10 > 4.30103)
sum(Day10Data$count_freq10 > 4.30103)
sum(Day14Data$count_freq10 > 4.30103)

mean(Day7Data$count_freq10)
mean(Day10Data$count_freq10)
mean(Day14Data$count_freq10)


