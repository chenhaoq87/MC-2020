#Checking on sensitivity analysis proportions
#by SL on 9/8/20

library(dplyr)
sensday7<- data3 %>% filter(day == 7)
sum(sensday7$count > 4.30103)

sensday10<- data3 %>% filter(day == 10)
sum(sensday10$count > 4.30103)

sensday14<- data3 %>% filter(day == 14)
sum(sensday14$count > 4.30103)