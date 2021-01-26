#initial microbial count

#Sam Lau 
#02/14/20

library(readr)
dat <- read_csv("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Initial Microbial Count/Analysis/Initial microbial count with filtering_rawdata_021420.csv")




#dat <- read_csv("/Volumes/ag-food/FOOD-Labs/BoorWiedmannLab/sl2763_Samantha Lau/Projects/2019 FFAR project/Initial Microbial Count/Analysis/dat_refined_062419.csv")

#Convert dates to appropriate format
date_processed<- as.Date(dat$dateProcessed,'%m/%d/%Y')
date_initial<- as.Date(dat$dateInitial,'%m/%d/%Y')

diff_in_days<- difftime(date_initial ,date_processed, units = "days") # days

#Make a column called Day_Initial_Actual
dat$Day_Initial_Actual<-paste(diff_in_days)

#Make a df only including ID, dateProcessed, dateInitial, Day_Initial_Actual, and all the DI, D7, D10, D14 data; removed "countability" columns
dat2 <- dat[c(1,21,9:10,12:13,15:16,18:19)]

#Make df without day 14
#dat3 <- dat2[c(1:8)]

#include day 14
dat3 <- dat2
#Convert characters to factors
dat3 <- as.data.frame(unclass(dat3))

df <- data.frame(dat3, dat3[6])
colnames(df)[colnames(df)=="D7datavalue.1"] <- "D7_data"

#If it says NA, then we put in the column the data value because that means the data value is correct
#Otherwise it would just be < or >
#Duplicate the values
#In order to get the parameters that go into the line, you want the actual data 
#If its not NA in that column, it will become NA in the value column

df$D7greaterthanlessthan<- as.character(df$D7greaterthanlessthan)
df$D7_data<- as.character(df$D7_data)
df$D7greaterthanlessthan[is.na(df$D7greaterthanlessthan)] <- df$D7_data[is.na(df$D7greaterthanlessthan)]

df$D7greaterthanlessthan[df$D7greaterthanlessthan == "<"] <- NA
df$D7greaterthanlessthan[df$D7greaterthanlessthan == ">"] <- NA
colnames(df)[colnames(df)=="D7greaterthanlessthan"] <- "D7"


##day 10
df <- data.frame(df, dat3[8])
colnames(df)[colnames(df)=="D10datavalue.1"] <- "D10_data"

df$D10greaterthanlessthan<- as.character(df$D10greaterthanlessthan)
df$D10_data<- as.character(df$D10_data)
df$D10greaterthanlessthan[is.na(df$D10greaterthanlessthan)] <- df$D10_data[is.na(df$D10greaterthanlessthan)]

df$D10greaterthanlessthan[df$D10greaterthanlessthan == "<"] <- NA
df$D10greaterthanlessthan[df$D10greaterthanlessthan == ">"] <- NA
colnames(df)[colnames(df)=="D10greaterthanlessthan"] <- "D10"


##day initial
df <- data.frame(df, dat3[4])
colnames(df)[colnames(df)=="DIdatavalue.1"] <- "DI_data"

df$DIgreaterthanlessthan<- as.character(df$DIgreaterthanlessthan)
df$DI_data<- as.character(df$DI_data)
df$DIgreaterthanlessthan[is.na(df$DIgreaterthanlessthan)] <- df$DI_data[is.na(df$DIgreaterthanlessthan)]

df$DIgreaterthanlessthan[df$DIgreaterthanlessthan == "<"] <- NA
df$DIgreaterthanlessthan[df$DIgreaterthanlessthan == ">"] <- NA
colnames(df)[colnames(df)=="DIgreaterthanlessthan"] <- "DI"

##day 14
df <- data.frame(df, dat3[10])
colnames(df)[colnames(df)=="D14datavalue.1"] <- "D14_data"

df$D14greaterthanlessthan<- as.character(df$D14greaterthanlessthan)
df$D14_data<- as.character(df$D14_data)
df$D14greaterthanlessthan[is.na(df$D14greaterthanlessthan)] <- df$D14_data[is.na(df$D14greaterthanlessthan)]

df$D14greaterthanlessthan[df$D14greaterthanlessthan == "<"] <- NA
df$D14greaterthanlessthan[df$D14greaterthanlessthan == ">"] <- NA
colnames(df)[colnames(df)=="D14greaterthanlessthan"] <- "D14"


##make new cleaner df by removing the data value columns
df2 <- df[c(1:3,5,7,9)]

#make a new row with number of NAs
df2$na_count <- apply(df2, 1, function(x) sum(is.na(x)))
df2$na_count <- as.factor(df2$na_count)

#want to see the number of NAs in the columns
#want the number of NAs to be less than or equal to 2
#so that you have at least two points that can be used to get equation of line
summary(df2)

#data frame without the NA count
#df3<- df2[c(1,3:7)]

#create a new df 
#have columns be ID, Day, and count
#need to flip the column names
#need the column ID to be listed as many times as there are Day values


#maybe create a df with day 0, 7, 10, 14 automatically put in as the days
#any of the rows with NA in that specific day, delete later
#maybe replicate the ID row so it appears 4 times
#want to expand the df such that each row appears the number of times specified in the column NA
#create new column with 4 subtract na_count to write the frequency that ID should appear
df3<-df2
df3['Four'] = '4'
df3$Four <- as.numeric(as.character((df3$Four)))


df3$na_count <- as.numeric(as.character((df3$na_count)))
#df3$na_count<-(df3$na_count-1)
df3$datapoints<- (df3$Four-df3$na_count)
df3$datapoints <- as.numeric(as.character((df3$datapoints)))

#create a new data frame without the na_count and the Four column
df4<-df3[c(1:6,9)]

#delete rows with <2 data points DIDNT WORK and make it a new dataframe
#df4[!df4$datapoints == "0"]
#subset(df4, datapoints==0)
df5 <- df4[df4$datapoints>=2,]

#need to create the day column with day initial, 7, 10, 14 for each ID
#for each ID, add rows with the value in the second column being 7 10 or 14
#maybe an if statement with adding rows? 
# example, if d7 column has a value (or not NA) then add a row with the same ID and 7
#or create a new df with all the d7, all the d10, and d14 and then bind them later

day7<-df5[c(1,4)]
day7['Day'] = '7'
day7<-day7[c(1,3,2)]
day7$D7 <- as.numeric(as.character((day7$D7)))
day7B <- day7[day7$D7>=0,]
day7revised <- na.omit(day7B)
colnames(day7revised)[colnames(day7revised)=="D7"] <- "Count"

#do the same with day10
day10<-df5[c(1,5)]
day10['Day'] = '10'
day10<-day10[c(1,3,2)]
day10$D10 <- as.numeric(as.character((day10$D10)))
day10B <- day10[day10$D10>=0,]
day10revised <- na.omit(day10B)
colnames(day10revised)[colnames(day10revised)=="D10"] <- "Count"

#do the same with day14
day14<-df5[c(1,6)]
day14['Day'] = '14'
day14<-day14[c(1,3,2)]
day14$D14 <- as.numeric(as.character((day14$D14)))
day14B <- day14[day14$D14>=0,]
day14revised <- na.omit(day14B)
colnames(day14revised)[colnames(day14revised)=="D14"] <- "Count"

#do the same with DI
dayI<-df5[c(1:3)]
dayI$DI <- as.numeric(as.character((dayI$DI)))
dayIB <- dayI[dayI$DI>=0,]
dayIrevised <- na.omit(dayIB)
#dayIrevised(Day_Initial_Actual)[dayIrevised(Day_Initial_Actual)=="Day_Initial_Actual"] <- "Day"
colnames(dayIrevised)[colnames(dayIrevised)=="Day_Initial_Actual"] <- "Day"
colnames(dayIrevised)[colnames(dayIrevised)=="DI"] <- "Count"

#bind all data frames
#Initialdatacombined<-as.data.frame(rbind(dayIrevised, day7revised, day10revised, day14revised))
Initialdatacombined<-(rbind(dayIrevised, day7revised, day10revised, day14revised))
#Initialdatacombined<-Initialdatacombined[c(1:3)]
sorted.out <- Initialdatacombined[order(as.numeric(as.character(Initialdatacombined$ID))), ]
#need to check if day initial is included in the spreadsheet
write.csv(sorted.out, "Look_at_sorted_data_SL_021420.CSV")

#Create code that returns FALSE or TRUE if it has day initial counts
#For loop
#In the beginning part of the code before it goes and make a line
#Make a for loop that pushes out either a table or dataframe with the unique ID and return true or false
#True if it has day initial
#False if it doesn't
#Then write code that counts up the amount of trues and the amount of falses

sampleID_list <- unique(sorted.out$ID)

#create a list to put the true/false results in, after all the loops we will bind them
#and write it to a file
truefalse_list <- list()

#we have to keep track of our list position so we can add the new row to the end
list_index <- 1 ##first element will go in 1, then we will increase so next goes in 2, etc
sampleID <- sampleID_list[[1]]
#loop through each ID
for (sampleID in sampleID_list) {

  this_id_only <- sorted.out[sorted.out$ID == sampleID, ]
  Day_values<-as.numeric(as.character(this_id_only$Day))

  TRUEFALSE <- 0

  if (min(Day_values) < 7){
    TRUEFALSE <- 1}
  
  
  #if (min(Day_values<7)){
   # TRUEFALSE<-0}
  #or do TRUEFALSE<-YES
  #else {TRUEFALSE<-1
  #or do TRUEFALSE<-NO
  
  
  
  
  
  #start making a new df with all the results
  #want a separate column with ID,slope and Y-intercept
  truefalseresults <- list (sampleID = sampleID, 
                   TRUEFALSE = TRUEFALSE)
  
  
  #now add these results to the end of our list of results
  truefalse_list[[list_index]] <- truefalseresults
  #remember we need to increase i by 1 so the next time we add our new results to the end
  list_index <- list_index + 1
  
} 


#bind all the items of our list together, 
#the do.call function applies rbind to each item of the list
#convert this to a dataframe for proper output
list1 <- data.frame(do.call(rbind, truefalse_list))

#write.csv(output, "Initial_microbial_contamination_SL_082919.csv")

#Make all the elemnts in the dataframe a single character?
list2 = data.frame(lapply(list1, as.character), stringsAsFactors=FALSE)

# write file
#write.csv(list2,"TrueFalse_SL_090519.csv")




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
Cleaneddata<-sorted.out
Cleaneddata$Count<-log10(Cleaneddata$Count)
ID_list <- unique(Cleaneddata$ID)

#create a list to put results in, after all the loops we will bind them
#and write it to a file
output_list <- list()

#we have to keep track of our list position so we can add the new row to the end
list_index <- 1 ##first element will go in 1, then we will increase so next goes in 2, etc
ID <- ID_list[[1]]
#loop through each ID
for (ID in ID_list) {
  
  this_id_only <- Cleaneddata[Cleaneddata$ID == ID, ]
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

# write file
write.csv(output2,"Initialmicrodata_final_SL_021420.csv")


# Now we just to write our list to a file
#THIS DIDN'T WORK
#write.table(output, file = "Variables_for_equation_SL_082919.csv", sep=",", row.names=FALSE)


