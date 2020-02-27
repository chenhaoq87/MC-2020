#Monte Carlo Model Validation using Kolmogorov-Smirnov test
library(dplyr)
library(ggplot2)
Predicted_counts <- subset(data,data$day == 14 & data$count >1)
VSL_import <- read.csv("VSL_Gpos.csv")
VSL_counts <- VSL_import$log_VSL_D14
ks.test(Predicted_counts$count, VSL_counts)
hist(VSL_counts, col=rgb(0.1, 0.1, 0.1, 0.5), xlim=c(1,8), 
     ylim=c(0,200), xlab = expression(paste(LOG[10],CFU/mL)))
hist(Predicted_counts$count, col=rgb(0.8, 0.8, 0.8, 0.5), add=T)

ggplot(Predicted_counts,aes(x=count,y=..density..)) +
  #geom_histogram(alpha=0.5) +
  #geom_histogram(data=data.frame(count=VSL_counts),
   #              aes(x=count,y=..density..),
  #               fill="red",
  #               alpha=0.5) +
  stat_density(alpha=0.5,n=512,
               bw=.75) +
  stat_density(data=data.frame(count=VSL_counts),
               fill="red",
               alpha="0.5",
               bw=.75,
               n=512)

data.frame(rbind(cbind(count=Predicted_counts$count,proj="monte"),
                                   cbind(count=VSL_counts,proj="VSL"))) -> for_boxplots

for_boxplots %>%
  mutate(count=as.numeric(as.character(count))) %>%
  filter(!is.na(count)) -> for_boxplot
#data.frame(rbind(cbind(count=Predicted_counts$count,proj="monte"),
                 #cbind(count=VSL_counts,proj="VSL"))) 



ggplot(for_boxplot,aes(x=proj,y=count)) + geom_boxplot(alpha=0.5) + 
  theme_bw() + coord_flip() + theme(text=element_text(size=14, family = "Times New Roman"))

ggplot(for_boxplot,aes(x=count)) + stat_ecdf(aes(linetype=proj)) + 
  theme_bw()+ theme(text=element_text(size=14, family = "Times New Roman"))+
  labs(x=expression(paste(LOG[10],CFU/mL)), y= "Density")

for_boxplot %>% filter(proj=="monte") %>% select(count) -> monte
for_boxplot %>% filter(proj=="VSL") %>% select(count) -> vsl
