#Amt of extra days

#check to see if any start at below 1 cfu/half gallon
data4<-data3
sum(data4$logMPN_init < -3.276921)

lessthan1<- subset(data4, data4$logMPN_init < -3.28)
lessthan1$count <- -8
greaterthan1<- subset(data4, data4$logMPN_init > -3.28)

#combine the 2 datasets
combined<-rbind(lessthan1, greaterthan1)

#checking proportions

library(dplyr)
sensday7<- combined %>% filter(day == 7)
sum(sensday7$count > 4.30103) /100000
sensday10<- combined %>% filter(day == 10)
sum(sensday10$count > 4.30103)/100000
mean(sensday7$count)
mean(sensday10$count)


sensday7<- combined %>% filter(day == 7)
sum(sensday7$count > 4.30103) /100000
sensday8<- combined %>% filter(day == 8)
sum(sensday8$count > 4.30103) /100000
sensday9<- combined %>% filter(day == 9)
sum(sensday9$count > 4.30103) /100000
sensday10<- combined %>% filter(day == 10)
sum(sensday10$count > 4.30103)/100000
sensday11<- combined %>% filter(day == 11)
sum(sensday11$count > 4.30103)/100000
sensday12<- combined %>% filter(day == 12)
sum(sensday12$count > 4.30103)/100000
sensday13<- combined %>% filter(day == 13)
sum(sensday13$count > 4.30103)/100000
sensday14<- combined %>% filter(day == 14)
sum(sensday14$count > 4.30103)/100000


sensday7<- combined %>% filter(day == 7)
sum(sensday7$count > 4.30103) /100000
sensday10<- combined %>% filter(day == 10)
sum(sensday10$count > 4.30103)/100000
mean(sensday7$count)
mean(sensday10$count)



library(dplyr)
sensday7<- data3 %>% filter(day == 7)
sum(sensday7$count > 4.30103)

sensday10<- data3 %>% filter(day == 10)
sum(sensday10$count > 4.30103)

sensday14<- data3 %>% filter(day == 14)
sum(sensday14$count > 4.30103)

mean(sensday7$count)
mean(sensday8$count)
mean(sensday9$count)
mean(sensday10$count)
mean(sensday11$count)
mean(sensday12$count)
mean(sensday13$count)
mean(sensday14$count)