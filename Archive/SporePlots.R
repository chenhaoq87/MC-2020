# Author: Mike Phillips
# Date: Thu Jan 16 10:03:56 2020
#  mdp38@cornell.edu

library(dplyr)
library(ggplot2)
library(rstudioapi)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext() 
setwd(dirname(current_path$path ))
print( getwd() )
source("MCsimV3_0.R")


#Analsis and plotting ----

data %>%
  group_by(day) %>%
  summarise(Mean = mean(count), SD = sd(count))-> mean_by_day

mean_by_day$day <- factor(mean_by_day$day)  #Requires discrete x values should be factors

#A ggplot Bar Chart
#THE BASICS OF HOW U DO A NICER BAR CHART
fig_ssc <- ggplot() + geom_bar(aes(y = Mean, x = day), data = mean_by_day, 
                               stat = "identity", color = "black", position = "dodge")

#Add the number above the bar plot                            
fig_ssc <- fig_ssc + geom_text(data = mean_by_day, 
                               aes(x = day, y = Mean, label = round(Mean, 2)), 
                               position=position_dodge(width = 1), 
                               vjust=-0.25, size = 4)

##Some optional stuff that shows how to customize the look of the plot
fig_ssc <- fig_ssc + scale_fill_grey(start = 0.8, end = 0.0) + theme_bw() +
  theme(legend.title = element_blank(), legend.position = "right",
        plot.title = element_text(face = "bold", size = 18,hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 14 ),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.ticks.x.bottom = element_blank())

## Set up the title of the graphs and the way the x and y axes are handled
fig_ssc <- fig_ssc + labs(y = "Mean Bacteria Count (log cfu/ml)", x = "Day" )+
  scale_x_discrete(breaks = seq(0, 24, 1)) +
  scale_y_continuous(limits = c(-3, 9), expand = c(0, 0)) +
  ggtitle("Mean Bacteria Per Day")    

#type the name to display it
fig_ssc


days <- seq(1, 24, 1)
spoiled <- vector(mode = "logical", length(days))
spoiled_by_day <- data.frame(days, spoiled)
for (d in days) {
  numerator <-   length(which(subset(data, data$day == d)$count>4.3))
  denominator <- length(subset(data, data$day == d)$count)
  
  spoiled_by_day[spoiled_by_day$days==d,]$spoiled <- numerator / denominator * 100
}

spoiled_by_day$days <- factor(spoiled_by_day$days)
fig_SBD <- ggplot() + geom_bar(aes(y = spoiled_by_day$spoiled, x = spoiled_by_day$days), data = spoiled_by_day, 
                               stat = "identity", color = "black", position = "dodge")


#Add the number above the bar plot                            
fig_SBD <- fig_SBD + geom_text(data = spoiled_by_day, 
                               aes(x = spoiled_by_day$days, y = spoiled_by_day$spoiled, label = round(spoiled_by_day$spoiled, 2)), 
                               position=position_dodge(width = 1), 
                               vjust=-0.25, size = 4)


##Some optional stuff that shows how to customize the look of the plot
fig_SBD <- fig_SBD + scale_fill_grey(start = 0.8, end = 0.0) + theme_bw() +
  theme(legend.title = element_blank(), legend.position = "right",
        plot.title = element_text(face = "bold", size = 18,hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 14 ),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.ticks.x.bottom = element_blank())


## Set up the title of the graphs and the way the x and y axes are handled
fig_SBD <- fig_SBD + labs(y = "Percentage > 4.3log", x = "Day" )+
  scale_x_discrete(breaks = seq(1, 24, 1), labels=seq(1, 24, 1)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  ggtitle("Percent of samples spoiled by day")    

#type the name to display it
fig_SBD
