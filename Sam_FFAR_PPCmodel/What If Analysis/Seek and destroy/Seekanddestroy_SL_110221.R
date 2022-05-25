load("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Initial Microbial Count/Analysis/CheckingInitialDistribution/What if analysis with 1.11log initial count/Seek and destroy/BASEMODEL_ppcmodel_1.11initial.RData")

#(b) Calculate proportion of samples that have spoiled WHEN Remove ST 13, the ST with the highest growth rate by setting final conc as -8
#+ The proportions of the other ST will remain unchanged

#Set the final concentration of ST 100 as -8
data6<-data3
data6$count <- ifelse(data6$AT=="100",-8,data6$count)

sensday7b<- data6 %>% filter(day == 7)
sensday8b<- data6 %>% filter(day == 8)
sensday9b<- data6 %>% filter(day == 9)
sensday10b<- data6 %>% filter(day == 10)
sensday14b<- data6 %>% filter(day == 14)

sum(sensday7b$count > 4.30103)/100000
sum(sensday8b$count > 4.30103)/100000
sum(sensday9b$count > 4.30103)/100000
sum(sensday10b$count > 4.30103)/100000
sum(sensday14b$count > 4.30103)/100000

mean(sensday7b$count)
mean(sensday8b$count)
mean(sensday10b$count)
mean(sensday14b$count)



#(d) Calculate proportion of samples when setting the final concentration of ST 13 and 9 as 0
#+ The proportions of the other ST will remain unchanged

#set final concentration of ST 13 and 9 as zero, the two most frequent ST
data7<-data3
data7$count <- ifelse(data7$AT=="13",-8,data7$count)
data7$count <- ifelse(data7$AT=="9",-8,data7$count)

sensday7c<- data7 %>% filter(day == 7)
sensday8c<- data7 %>% filter(day == 8)
sensday9c<- data7 %>% filter(day == 9)
sensday10c<- data7 %>% filter(day == 10)
sensday14c<- data7 %>% filter(day == 14)

sum(sensday7c$count > 4.30103)/100000
sum(sensday8c$count > 4.30103)/100000
sum(sensday9c$count > 4.30103)/100000
sum(sensday10c$count > 4.30103)/100000
sum(sensday14c$count > 4.30103)/100000

mean(sensday7c$count)
mean(sensday9c$count)
mean(sensday10c$count)
mean(sensday14c$count)

