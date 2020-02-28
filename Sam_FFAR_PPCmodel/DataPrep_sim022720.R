library(readr)
library(tidyverse)

# load data
library(readxl)
df<- read_excel("Sam_FFAR_PPCmodel/Samantha Sample micro data 2020-02-25.xlsx")
freq_df <- read_csv("Sam_FFAR_PPCmodel/InputFiles/Frequency_ALLISOLATES_021120.csv") #freq file; isolates from Sam R's study
# df <- read_csv("Sam_FFAR_PPCmodel/Samantha Sample micro data 2020-02-25.csv") #count data for samples from Sam R's study

# rename columns in df
colnames(df)[1] <- "VSLNumber"
# only use pasteurized fluid white milk; either 1, 2, 3, or 4 from the right (after removing "-" to the left)
df$SampleID<- trimws(df$SampleID, which = c("both"))
df$temp <- df$SampleID
df <- df %>%
  separate(temp,c("id1","id2"),sep="-")
df <- df %>%
  separate(id1,c("id3","milkType"),sep=-2)
dfw <- df %>%
  filter(milkType == "1"|milkType == "2"|milkType == "3"|milkType == "4") %>% #only include white milk
  filter(id3 !="1N 4"&id3 !="1N 3"&id3 !="1N 2"&id3 !="1N 1"&id3 !="2N 1"&id3 !="3NPOST 1"&id3 !="3NPRE 1" )%>%
  filter(SampleID!=	"72-1*" & SampleID!=	"73-1*")

# prep freq_df
# freq_df$VSLNumber_SampleID <- as.factor(freq_df$VSLNumber_SampleID)
# summary(freq_df$VSLNumber_SampleID)
# freq_df1 <- freq_df
# # get milkType and ids from freq_df
# freq_df$temp <- freq_df$VSLNumber_SampleID
# freq_df <- freq_df %>%
#   separate(temp,c("VSLNumber","SampleID"),sep="_")
# freq_df$VSLNumber <- as.numeric(freq_df$VSLNumber)
# freq_df$SampleID2 <- freq_df$SampleID
# freq_df <- freq_df %>%
#   separate(SampleID2,c("id3","id4"),sep="-")
# freq_df <- freq_df %>%
#   separate(id3,c("id5","milkType"),sep=-1)
# freq_dfw <- freq_df %>%
#   filter(milkType == "1"|milkType == "2"|milkType == "3"|milkType == "4") #only include white milk; result is same since all in freq_df were white milk

# # review df Sam R had sent to make sure OK
freq_df_check <- read_csv("Sam_FFAR_PPCmodel/Copy of 16S ST by Sample for Sarah 2020-02-21.csv")
freq_df_check$VSLNumber <- freq_df_check$`VSL Number`
freq_df_check$VSLNumber<-as.numeric(freq_df_check$VSLNumber)
freq_df_check$temp <- freq_df_check$SampleID
freq_df_check <- freq_df_check %>%
  separate(temp,c("id1","id2"),sep="-")
freq_df_check <- freq_df_check %>%
  separate(id1,c("id3","milkType"),sep=-2)
freq_df_checkw <- freq_df_check %>%
  filter(milkType == "1"|milkType == "2"|milkType == "3"|milkType == "4") #only include white milk
colnames(freq_df_checkw)[5] <- "X16S_ST"
freq_df_checkw <- freq_df_checkw %>%
  filter(id3 !="1N 4"&id3 !="1N 3"&id3 !="1N 2"&id3 !="1N 1"&id3 !="2N 1"&id3 !="3NPOST 1"&id3 !="3NPRE 1" )%>%
  filter(SampleID!=	"72-1*" & SampleID!=	"73-1*")
freq_df_checkw$VSLNumber_SampleID <- paste(freq_df_checkw$VSLNumber,freq_df_checkw$SampleID,sep="_")
freq_df_checkw_unique <- unique(freq_df_checkw[c(1,6,3,5,10)])
# set1 <- freq_df_checkw_unique[c(2:4)]
# set2 <- freq_df[c(6,7,5)]
# temp <- setdiff(set1,set2) # why missing 16S = 42 for sample 22-2 in VSL project 457
# temp2 <- setdiff(set2,set1) #ok

# there are a number of instances that have repeat STs for 1 sample; not unique?
# set2_sum <- set2 %>% 
#   group_by(X16S_ST,VSLNumber,SampleID)%>%
#   summarize(n=n())

# Repeat STs were likely isolates from different days in shelf life and had different Genus_species & were Pseudomonas; YES; use freq_df without changing

#clear environment
rm(list=setdiff(ls(), c("dfw","freq_df1")))

# prep dfw
dfw$VSLNumber_SampleID <- paste(dfw$VSLNumber,dfw$SampleID,sep="_")
dfw$Test_Day <- paste(dfw$Test,dfw$Day,sep="_")
dfw$Plant_VSLNumber_SampleID <- paste(dfw$Plant,dfw$VSLNumber_SampleID,sep="_")
dfw$Countability[is.na(dfw$Countability)==TRUE] <- "OK"
dfw <- dfw %>%
  filter(Countability!="M")
# convert dates to appropriate format
dfw$`Date Collected` <- as.Date(dfw$`Date Collected`,'%m/%d/%Y')
dfw$`Date Processed`<- as.Date(dfw$`Date Processed`,'%m/%d/%Y')
dfw$Day_Initial_Actual<- difftime(dfw$`Date Collected`,dfw$`Date Processed`, units = "days") # days

#temp
# dfw$log10count<- log10(as.numeric(dfw$`Count Value`))
colnames(dfw)[9]<-"greaterThanLessThan"
# dfw2 <- dfw[c(16,19,17)] %>%
#   spread(.,Test_Day,log10count)"
dfw2 <- dfw[c(16,10,17)] %>%
  spread(.,Test_Day,`Count Value`)
dfw2[is.na(dfw2)==TRUE] <- "NotTested"

dfw$LODtemp <- "LOD"
dfw$LOD_Test_Day <- paste(dfw$LODtemp,dfw$Test_Day,sep="_")
dfw2b <- dfw[c(17,9,20)] %>%
  spread(.,LOD_Test_Day,greaterThanLessThan)
dfw2b[is.na(dfw2b)==TRUE] <- "OK"
dfw2_merged <- merge(dfw2,dfw2b,by="Plant_VSLNumber_SampleID")
temp <- unique(dfw[c(17,18)])
dfw2_merged2 <- merge(dfw2_merged,temp,by="Plant_VSLNumber_SampleID")
dfw2_merged2<- separate(data=dfw2_merged2, Plant_VSLNumber_SampleID,c("Plant", "VSLNumber", "SampleID"),sep="_")
dfw2_merged2 <- dfw2_merged2[c(1:3,28,4:27)]
# write.csv(dfw2_merged2,"dfw2_merged_022720.csv")

# read in Sam R's PPC data
df_PPC_raw <- read_excel("Sam_FFAR_PPCmodel/Sarah S-PPC query 2020-02-27_withPlant.xlsx")
# rename columns in df_PPC
df_PPC <- df_PPC_raw  
colnames(df_PPC)[2] <- "VSLNumber"
# only use pasteurized fluid white milk; either 1, 2, 3, or 4 from the right (after removing "-" to the left)
df_PPC$SampleID<- trimws(df_PPC$SampleID, which = c("both"))
df_PPC$temp <- df_PPC$SampleID
df_PPC <- df_PPC %>%
  separate(temp,c("id1","id2"),sep="-")
df_PPC <- df_PPC %>%
  separate(id1,c("id3","milkType"),sep=-2)
df_PPCw <- df_PPC %>%
  filter(milkType == "1"|milkType == "2"|milkType == "3"|milkType == "4") %>% #only include white milk
  filter(id3 !="1N 4"&id3 !="1N 3"&id3 !="1N 2"&id3 !="1N 1"&id3 !="2N 1"&id3 !="3NPOST 1"&id3 !="3NPRE 1" )%>%
  filter(SampleID!=	"72-1*" & SampleID!=	"73-1*")
temp2 <- df_PPCw[c(1:6)]
# merge Sam R PPC with dfw2merged2
SamRdat_full <- merge(dfw2_merged2,temp2,by=c("Plant","VSLNumber","SampleID"), all.x=TRUE)
# write.csv(SamRdat_full,"SamRdat_full.csv")
SamRdat_full$`S-PPC` <- as.character(SamRdat_full$`S-PPC`)
#clear environment
# rm(list=setdiff(ls(), c("SamRdat_full","freq_df1")))
test <- unique(freq_df[c(1,5)]) # unique ST
df_freq <- read_excel("Sam_FFAR_PPCmodel/Samantha Sarah Frequency Table 2020-02-26.xlsx")
df_freq2 <- df_freq %>%
  filter(Day == "D21") %>%
  filter(`Gram Status`== "-")
df_freq2$temp <- df_freq2$SampleID
df_freq2 <- df_freq2 %>%
  separate(temp,c("id1","id2"),sep="-")
df_freq2 <- df_freq2 %>%
  separate(id1,c("id3","milkType"),sep=-2)
df_freq2w <- df_freq2 %>%
  filter(milkType == "1"|milkType == "2"|milkType == "3"|milkType == "4") %>% #only include white milk
  filter(id3 !="1N 4"&id3 !="1N 3"&id3 !="1N 2"&id3 !="1N 1"&id3 !="2N 1"&id3 !="3NPOST 1"&id3 !="3NPRE 1" )%>%
  filter(SampleID!=	"72-1*" & SampleID!=	"73-1*")
df_freq2w2 <- unique(df_freq2w[c(4,9,10,16)])
df_freq2w3 <- df_freq2w2 %>%
  filter(is.na(`16S ST`)!=TRUE)
colnames(df_freq2w3)[1] <- "Plant"

#
# write.csv(df_freq2w3,"df_freq_D21only_unique16SST.csv")

## get data for Upstate Buffalo
UB_16Sfreq <- df_freq2w3 %>%
  filter(Plant == "Upstate Buffalo")
UB_micro_all <- SamRdat_full %>%
  filter(Plant == "Upstate Buffalo")
UB_micro_S_PPC <- UB_micro_all %>%
  filter(`S-PPC` == "TRUE")
UB_micro_S_PPC <- UB_micro_S_PPC[c(1:4,10,9,5:8,16,15,11:14,17:31)]
UB_micro_all <- UB_micro_all[c(1:4,31,10,9,5:8,16,15,11:14,17:30)]

# 
