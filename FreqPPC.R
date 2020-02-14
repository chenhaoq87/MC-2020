#load packages
library(readr)
library(dplyr)
library(tidyr)

## prepare VSL data for frequency of PPC; here defined as whether or not sample exceed 20,000 cfu/mL SPC by d10 (called "fail")
vsl_2015aTo2017b <- read_csv("sim3_VSL_020320.csv")
vsl_2015aTo2017b <- vsl_2015aTo2017b %>%
  filter(dataValue!="NT") #remove rows where "NT" (i.e., not tested)
vsl_2015aTo2017b$log10SPC <- log10(as.numeric(vsl_2015aTo2017b$dataValue))

temp_spc <- vsl_2015aTo2017b[c(3,13,19)] %>%
  spread(.,Day,log10SPC) 

df_freqPPC <- temp_spc[c(1)]
df_freqPPC$exceed20000_d10 <- ifelse(temp_spc$DI>4.3 | temp_spc$D7>4.3 | temp_spc$D10>4.3, 1, 0)
df_freqPPC$exceed20000_d10[is.na(df_freqPPC$exceed20000_d10)] <- 0 # replace NA with zero
df_freqPPC <- df_freqPPC%>%
  separate(ID, c("plantID", "sampling","sampleID"), "_", remove = TRUE)
df_freqPPC <- as.data.frame(unclass(df_freqPPC))
df_freqPPC$exceed20000_d10 <- as.factor(df_freqPPC$exceed20000_d10)

#clean up environment
rm(vsl_2015aTo2017b,temp_spc)
