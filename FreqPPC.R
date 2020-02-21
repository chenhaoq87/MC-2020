#load packages
library(readr)
library(dplyr)
library(tidyr)

## prepare VSL data for frequency of PPC; here defined as whether or not sample exceed 20,000 cfu/mL SPC by d10 (called "fail")
vsl_2015aTo2017b <- read_csv("sim3_VSL_020320.csv")
vsl_2015aTo2017b <- vsl_2015aTo2017b %>%
  filter(dataValue!="NT") #remove rows where "NT" (i.e., not tested)

vsl_2015aTo2017b$log10SPC <- log10(as.numeric(vsl_2015aTo2017b$dataValue))

#replace data below LOD with 25% of detection limit (10 CFU/mL so 2.5 CFU/mL)
vsl_2015aTo2017b$temp <- (vsl_2015aTo2017b$log10SPC)*0.25
vsl_2015aTo2017b$count_type[is.na(vsl_2015aTo2017b$count_type) == TRUE] <- "OK"
vsl_2015aTo2017b$log10SPC <- ifelse(vsl_2015aTo2017b$count_type == "belowLOD", vsl_2015aTo2017b$temp, vsl_2015aTo2017b$log10SPC)

#for data above the detection limit, for this purpose, exact count is not needed since above minimum of 4.3 log CFU/mL
temp_spc <- vsl_2015aTo2017b[c(3,13,19)] %>%
  spread(.,Day,log10SPC) 

df_freqPPC <- temp_spc[c(1)]
df_freqPPC$exceed20000_d10 <- ifelse(temp_spc$DI>4.3 | temp_spc$D7>4.3 | temp_spc$D10>4.3, 1, 0)
df_freqPPC$exceed20000_d10[is.na(df_freqPPC$exceed20000_d10)] <- 0 # replace NA with zero
df_freqPPC <- df_freqPPC%>%
  separate(ID, c("plantID", "sampling","sampleID"), "_", remove = TRUE)
df_freqPPC <- as.data.frame(unclass(df_freqPPC))
df_freqPPC$exceed20000_d10 <- as.factor(df_freqPPC$exceed20000_d10)

# ## Check how to do PPC stuff for Sam
# library(readxl)
# df <- read_excel("Sam_FFAR_PPCmodel/Initial microbial count with CVTA and SPC_112519.xlsx")
# # df2 <- df %>%
# #   filter(plantID=="36-1661") #plantI
# df <- df[c(1,9:32)]
# df$CVTA_DI <- ifelse(is.na(df$DIcvtagreaterthanlessthan == TRUE), df$DIcvtadatavalue, df$DIcvtagreaterthanlessthan)
# df$CVTA_D7 <- ifelse(is.na(df$D7cvtagreaterthanlessthan == TRUE), df$D7cvtadatavalue, df$D7cvtagreaterthanlessthan)
# df$CVTA_D10 <- ifelse(is.na(df$D10cvtagreaterthanlessthan == TRUE), df$D10cvtadatavalue, df$D10cvtagreaterthanlessthan)
# df$CVTA_D14 <- ifelse(is.na(df$D14cvtagreaterthanlessthan == TRUE), df$D14cvtadatavalue, df$D14cvtagreaterthanlessthan)
# df$CVTA_DI[df$DIcvtagreaterthanlessthan == ">"] <- 1
# df$CVTA_D7[df$D7cvtagreaterthanlessthan == ">"] <- 1
# df$CVTA_D10[df$D10cvtagreaterthanlessthan == ">"] <- 1
# df$CVTA_D14[df$D14cvtagreaterthanlessthan == ">"] <- 1
# df$SPC_DI <- ifelse(is.na(df$DIspcgreater == TRUE), df$DIspcdatavalue, df$DIspcgreater)
# df$SPC_D7 <- ifelse(is.na(df$D7spcgreater == TRUE), df$D7spcdatavalue, df$D7spcgreater)
# df$SPC_D10 <- ifelse(is.na(df$D10spcgreater == TRUE), df$D10spcdatavalue, df$D10spcgreater)
# df$SPC_D14 <- ifelse(is.na(df$D14spcgreater == TRUE), df$D14spcdatavalue, df$D14spcgreater)
# 
# df$log10SPC_DI <- log10(as.numeric(df$SPC_DI))
# df$log10SPC_D7 <- log10(as.numeric(df$SPC_D7))
# df$log10SPC_D10 <- log10(as.numeric(df$SPC_D10))
# df$log10SPC_D14 <- log10(as.numeric(df$SPC_D14))
# 
# # df2 <- df[c(1,26:33)]
# df$exceed20000_d10 <- ifelse(df$log10SPC_DI>4.3 | df$log10SPC_D7>4.3 | df$log10SPC_D10>4.3, 1, 0)
# df$exceed20000_d10 [is.na(df$exceed20000_d10)] <- 0 # replace NA with zero
# 
# df$log10CVTA_DI <- log10(as.numeric(df$CVTA_DI))
# df$log10CVTA_D7 <- log10(as.numeric(df$CVTA_D7))
# df$log10CVTA_D10 <- log10(as.numeric(df$CVTA_D10))
# df$log10CVTA_D14 <- log10(as.numeric(df$CVTA_D14))
# 
# df <- df %>%
#   filter(!is.na(log10SPC_D14))
# 
# df$GN <- ifelse(is.na(df$log10CVTA_DI)==TRUE & is.na(df$log10CVTA_D7)==TRUE & is.na(df$log10CVTA_D10)==TRUE & is.na(df$log10CVTA_D10)==TRUE, 0, 1)
# 
# df$combo <- df$GN + df$exceed20000_d10
# 
# ggplot(df,aes(x=combo,fill=as.factor(GN)))+
#   geom_bar()
#clean up environment
rm(vsl_2015aTo2017b,temp_spc)

