#initial microbial count

#Sam Lau 
#02/27/20

library(readr)
library(plyr)
library(MASS)

dat <- read_csv("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Data prep for validation/Redo of initial concentration/Samdataset_readyforinitialconc_022720.csv")




#want to have 3 columns with unique ID, day, and count
#create a new df with all the d7, all the d10, and d14 and then bind them later
dat2<-dat

#creating df with only day7 data
day7<-dat2[c(1,5)]
day7["Day"]="7"
day7<-day7[c(1,3,2)]
day7$CVTA_D7 <- as.numeric(as.character((day7$CVTA_D7)))
day7<- na.omit(day7)
colnames(day7)[colnames(day7)=="CVTA_D7"] <- "Count"

#creating df with only day10 data
day10<-dat2[c(1,6)]
day10["Day"]="10"
day10<-day10[c(1,3,2)]
day10$CVTA_D10 <- as.numeric(as.character((day10$CVTA_D10)))
day10<- na.omit(day10)
colnames(day10)[colnames(day10)=="CVTA_D10"] <- "Count"

#creating df with only day14 data
day14<-dat2[c(1,7)]
day14["Day"]="14"
day14<-day14[c(1,3,2)]
day14$CVTA_D14 <- as.numeric(as.character((day14$CVTA_D14)))
day14<- na.omit(day14)
colnames(day14)[colnames(day14)=="CVTA_D14"] <- "Count"

#creating df with only dayI data
dayI<-dat2[c(1,3,4)]
dayI$CVTA_DI <- as.numeric(as.character((dayI$CVTA_DI)))
dayI <- na.omit(dayI)
colnames(dayI)[colnames(dayI)=="Day_Initial_Actual"] <- "Day"
colnames(dayI)[colnames(dayI)=="CVTA_DI"] <- "Count"


#bind all data frames
Initialdatacombined<-(rbind(dayI, day7, day10, day14))
sorted.out <- Initialdatacombined[order(as.numeric(as.character(Initialdatacombined$vsl_sample))), ]
write.csv(sorted.out, "Look_at_sorted_data_SL_022720.CSV")




#Now want to make a new df with ID, slope, and intercept
#in order to do this, you need to create an equation
#thus, need to create a df with all the variables necessary in the equation
#for line of best fit, it can be calculated by averaging all the x values
#averaging all the y values
#subtracting the individual (x,y) points from the means
#the y values would be counts
#the x values would be days
#have a column with ID, x, y, mean x, mean y, x-mean, y-mean, (x-mean)*(y-mean), (x-mean)^2
#slope is going to be (x-mean)*(y-mean)/ (x-mean)^2
#y intercept will be meanY-(slope*meanX)

#first grab the data frame that I made with list of IDs, day, and counts
#want to log transform the y values, so the counts
#if you don't log transform the values, it wouldn't be a linear equation

#you have to remove all the 0s because taking the log of that would be infinity
Cleaneddata<-sorted.out
Cleaneddata$Count[Cleaneddata$Count==0] <-0.001
Cleaneddata$Count<-log10(Cleaneddata$Count)
ID_list <- unique(Cleaneddata$vsl_sample)

#create a list to put results in, after all the loops we will bind them
#and write it to a file
output_list <- list()

#we have to keep track of our list position so we can add the new row to the end
list_index <- 1 ##first element will go in 1, then we will increase so next goes in 2, etc
ID <- ID_list[[1]]
#loop through each ID
for (ID in ID_list) {
  
  this_id_only <- Cleaneddata[Cleaneddata$vsl_sample == ID, ]
  X_values<-as.numeric(as.character(this_id_only$Day))
  Y_values<-as.numeric(as.character(this_id_only$Count))
  MeanX<-mean(as.numeric(as.character((X_values))))
  MeanY<-mean(as.numeric(as.character((Y_values))))
  
  #slope is going to be (x-mean)*(y-mean)/ (x-mean)^2
  #calculate the DENOMINATOR of slope equation:
  #when you do an operation on a list of things, it will do it to each item in the list
  #all the X values will be subtracted by the mean X
  Subtract_xmean <- X_values - MeanX
  Square_xmean <- Subtract_xmean^2
  Sum_square_xmean<-sum(Square_xmean)
  
  #slope is going to be (x-mean)*(y-mean)/ (x-mean)^2
  #calculate the NUMERATOR of slope equation:
  #when you do an operation on a list of things, it will do it to each item in the list
  #want the sum of (x-mean)*(y-mean)
  #(x-mean) was done earlier
  #only need to do the (y-mean)
  #then multiply them together and sum it
  
  Subtract_ymean<- Y_values- MeanY
  Multiply_xymeans<- Subtract_xmean * Subtract_ymean
  Sum_multiply_xymeans<- sum(Multiply_xymeans)
  
  #now that you have numerator and denominator
  #Put it all together for the slope equation
  #slope is going to be (x-mean)*(y-mean)/ (x-mean)^2
  Slope<-Sum_multiply_xymeans / Sum_square_xmean
  
  
  #now calculate the Y-intercept
  #y intercept will be meanY-(slope*meanX)
  Y_intercept<-MeanY-(Slope * MeanX)
  
  
  #start making a new df with all the results
  #want a separate column with ID,slope and Y-intercept
  results <- list (ID = ID, 
                   Slope = Slope,
                   Y_intercept = Y_intercept)

  
  #now add these results to the end of our list of results
  output_list[[list_index]] <- results
  #remember we need to increase i by 1 so the next time we add our new results to the end
  list_index <- list_index + 1
  
} 
  #closes inner isolate loop

#bind all the items of our list together, 
#the do.call function applies rbind to each item of the list
#convert this to a dataframe for proper output
output <- data.frame(do.call(rbind, output_list))

#write.csv(output, "Initial_microbial_contamination_SL_082919.csv")

#Make all the elemnts in the dataframe a single character?
output2 = data.frame(lapply(output, as.character), stringsAsFactors=FALSE)


#now use that file and add the slope and y-intercept columns and add them to get the LOG concentration
#do 10^log concentration to convert it into CFU/mL

output2$Slope <- as.numeric(as.character(output2$Slope))
output2$Y_intercept <- as.numeric(as.character(output2$Y_intercept))
output2["LOG10 Initial"]<- output2$Slope + output2$Y_intercept
output2["Concentration (CFU/mL)"]<- 10^(output2$`LOG10 Initial`)

# write file
write.csv(output2,"Initialmicrodata_final_SL_022720.csv")

#distribution of initial micro count


#use the initial microbial count column
loginitialcount<- output2
loginitialcount<-loginitialcount$`LOG10 Initial`


#attempt to plot out initial microbial count data
#histogram with frequencies
hist(loginitialcount, 10, col="black", main="Initial microbial count histogram", 
     xlab="LOG CFU/mL", xlim=c(-4,5), ylim = c(0,10), freq = TRUE, breaks = 20)

#histogram with proportions
hist(loginitialcount, 10, col="black", main="Initial microbial count histogram", 
     xlab="LOG CFU/mL", xlim=c(-4,5), ylim = c(0,0.5), freq = FALSE)

plot(density(loginitialcount))

#test for normality just to make sure the data is actually normal
shapiro.test(loginitialcount)

#produce the histogram for the normally distributed data (normal)
#is my data just not normally distributed?
hist(loginitialcount,probability=T, main="Initial microbial count histogram",
     xlab="LOG10 CFU/mL", ylim = c(0,0.6), xlim=c(-4,5), breaks = 20)
#add a density curve
#lines(density(loginitialcount),col=2)

#create a normal curve
meandata<-mean(loginitialcount)
print(meandata)
stddata<-sqrt(var(loginitialcount))
print(stddata)
curve(dnorm(x, mean=meandata, sd=stddata), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")
meandata
stddata

