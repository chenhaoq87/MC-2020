## ------------------------------------------------------------------------
##################################################
#  Monte Carlo Simulation 
##################################################
## Project: PPC Base Model
## Script purpose: Simulate PPC for half gallon samples of milk 
##                 over a given number of days. This is the validation file
##################################################
#Edited by SL on 11-01-21
##################################################


## ------------------------------------------------------------------------
library(splitstackshape)
library(dplyr)
library(readr)
library(truncnorm)
library(rmutil)

library(tidyverse)
library(censReg)
library(fitdistrplus)
library(rmutil)
#check which packages are in your file
#list.functions.in.file("Z:\\sl2763_Samantha Lau\\Projects\\2019 FFAR project\\Data prep for validation\\Model\\MCsimV2_Samdata_Samonlyplants_VALIDATION_SL_061220.R", alphabetic = TRUE)

sink(file = "TempNeg60filecatstatements_1.11init.txt")


cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S",'\n')," package import",'\n')


## ----knitr_options , include=FALSE---------------------------------------
set.seed(42)


## ------------------------------------------------------------------------
muAtNewTemp <- function(newTemp, oldMu, oldTemp = 6, T0 = -4.15) {
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  
  return(newMu)
}


## ------------------------------------------------------------------------
lagAtNewTemp <- function (newTemp, oldLag, oldTemp = 6, T0 = -4.15) {
  numerator <- oldTemp -T0
  denom <- newTemp - T0
  newLag <- ( (numerator / denom)^2) * oldLag
  return(newLag)
}



## ------------------------------------------------------------------------
buchanan_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax){
  ans <- LOG10N0 + (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) *     log(10)/mumax)) * mumax * (t - lag)/log(10) + (t >= lag) * (t > (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * (LOG10Nmax - LOG10N0)
  return(ans)
}

gompertz_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10N0 + (LOG10Nmax - LOG10N0) * exp(-exp(mumax * exp(1) * 
                                                      (lag - t)/((LOG10Nmax - LOG10N0) * log(10)) + 1))
  return(ans)
}


baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax -LOG10N0)))
  return(ans)
}



## ------------------------------------------------------------------------
#wrapper function because it calls the proper version
log10N <- function(t, lag, mumax, LOG10N0, LOG10Nmax,model_name="buchanan") {
  if (model_name == "buchanan") {
    return(buchanan_log10N(t, lag, mumax, LOG10N0, LOG10Nmax) )
  }
  else if(model_name == 'baranyi') {
    return(baranyi_log10N(t, lag, mumax, LOG10N0, LOG10Nmax) )
  }
  else if(model_name == 'gompertz') {
    return(gompertz_log10N(t, lag, mumax, LOG10N0, LOG10Nmax) )
  }
  else {
    stop(paste0(model_name, " is not a valid model name. Must be one of buchanan, baranyi, gompertz"))
  }
}


## ------------------------------------------------------------------------
#input files
frequency_file <- "Frequency_ALLISOLATES_021120.csv"
growth_file <- "ST_percentidentity_MCinput_SL_032021.csv"
init_file <- "Initialmicrodata_SamplantsONLY_MCinput_SL_022720.csv"

cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S",'\n'), " input files read in",'\n')

## ------------------------------------------------------------------------

freq_import <- read.csv(frequency_file, stringsAsFactors = FALSE, header = TRUE)
freq_data = freq_import$X16S_ST
freq_vec<-as.vector(freq_data) #before running sensitivity analysis (below) put a hashtag in front of this



## ------------------------------------------------------------------------
#Size is for n_sim bulk tanks, n_half_gal half gallon lots, n_day days (14-24)
n_sim <-10000    #1000 is for testing and exploring, experiments require at least 10k
n_halfgal <-10
n_day <- 14
start_day <- 1

#Repeat each element of the sequence 1..n_sim.Bulk tank data (MC runs)
BT <- rep(seq(1, n_sim), each = n_halfgal * n_day)
#Repeat the whole sequences times # of times
half_gal <- rep(seq(1, n_halfgal), times = n_day * n_sim)
#Vector of FALSE
AT <- vector(mode="logical", n_sim * n_halfgal * n_day)
#Repeat the days for each simulation run
day <- rep(rep(seq(start_day, start_day+n_day-1), each = n_halfgal), times = n_sim)
count <- vector(mode = "logical", n_sim * n_halfgal * n_day)

#matrix with columns:
#  BT   half_gal    AT    day   count
data <- data.frame(BT, half_gal, AT, day, count)

cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S",'\n')," matrix", '\n')

## ------------------------------------------------------------------------
#Import growth parameter data
growth_import <-read.csv(growth_file, stringsAsFactors = FALSE)
colnames(growth_import)[1]<-c("ST")
colnames(growth_import)[6]<-c("model_name")
growth_import$model_name <- as.factor(growth_import$model_name)
growth_import$lag <- as.numeric(growth_import$lag)
growth_import$mu <- as.numeric(growth_import$mu)
growth_import$N0 <- as.numeric(growth_import$N0)
growth_import$Nmax <- as.numeric(growth_import$Nmax)

#Import initial count logMPN data
initialcount_import <- read.csv(init_file, stringsAsFactors = FALSE)
#MPN Column
initialcount_data = initialcount_import[,3]
#LOG MPN Column
initialcountlog_data = initialcount_import[,2]

cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S",'\n')," data import done", '\n')

## ------------------------------------------------------------------------
### i. Temperature distribution 
#fixed to assign to each unique half gallon; Sam to double check
temps <- vector()
temps <- rlaplace(n=(n_sim*n_halfgal),m=1.624,s=2.31) # will be using temps vector later
# looked at documentation, here: https://www.rdocumentation.org/packages/ExtDist/versions/0.6-3/topics/Laplace
#https://doi.org/10.4315/0362-028X-73.2.312 (source)
# for (i in 1:n_sim*n_halfgal){temps <- sample(repvector,replace = F)} #to delete when sure

temps <- rep(NA, n_sim*n_halfgal)
for (i in 1:(n_sim*n_halfgal)){
  number <- rlaplace(1,m=1.624,s=2.31)
  while (number > 15 | number < -1) {
    number <- rlaplace(1,m=1.624,s=2.31)
  }
  temps[i] <- number
}



cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S",'\n')," temp distribution done",'\n')

### ii. Initial contamination
logMPN_mean <- c(0.3817872) #day initial
logMPN_sd <- c(1.108859) #day initial
logMPN_samp = rtruncnorm(n_sim, b=3, mean=logMPN_mean, sd=logMPN_sd)
#logMPN_samp = rnorm(n_sim, logMPN_mean, logMPN_sd)
MPN_samp = 10^logMPN_samp
MPN_samp_halfgal = MPN_samp * 1900 #MPN per half gallon (1892.71 mL in half gallon)

cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S",'\n')," initial contamination done",'\n')


### iii. Sample distributions
MPN_init<-vector()
allele <- vector()
for (i in 1:n_sim){
  MPN_init_samp <-rep(rpois(n_halfgal, MPN_samp_halfgal[i]), times = n_day)
  MPN_init<-c(MPN_init, MPN_init_samp)
  allele_samp <- rep(sample(freq_vec, n_halfgal, replace = T), times = n_day) 
  allele <- c(allele, allele_samp)
}

#Convert MPN_init from half-gallon to mLs
MPN_init_mL <- MPN_init / 1900
#remove 0's from the data and replace with detection limit
MPN_init_mL[MPN_init_mL == 0] <- 0.01;

data$logMPN_init <- log10(MPN_init_mL) #Add initial logMPN to data frame
data$AT<-allele #Add in AT data

cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S",'\n')," sample distributions drawn",'\n')

## ------------------------------------------------------------------------
# to allow for temperature assignment to unique BT & half gal, make a new variable BT_hg
data$BT_hg <- paste(data$BT,data$half_gal,sep="_")

# assign temps to unique BT_hg 
df_temp <- data[7]
df_temp %>% unique() -> df_temp 
df_temp$newT <- temps  #unhashtag this if you AREN'T doing validation



data2 <- merge(df_temp,data,by="BT_hg")
####

#make a temporary df with all except day
df_temp2 <- data2[c(1:5)]
df_temp2 %>%
  unique() -> df_temp2

# just determine newLag and newMu here
for (i in 1:(n_sim *n_halfgal)){
  #Find row in growth parameter data that corresponds to allele sample
  allele_index <- which(growth_import$ST == df_temp2$AT[i]) 
  
  #calculate the new growth parameters using the square root model and our
  #sampled temperature
  newLag <- lagAtNewTemp(df_temp2$newT[i], growth_import$lag[allele_index])
  newMu <-  muAtNewTemp(df_temp2$newT[i], growth_import$mu[allele_index])

  newMu<-newMu*0.684
  
  df_temp2$newLag[i] <- newLag
  df_temp2$newMu[i] <- newMu

}

####
# make a temporary df for merging
df_temp3 <- data[c(4:7)]

# now merge 
data3 <- merge(df_temp2,df_temp3,by="BT_hg")

#now determine count
for (i in 1:(n_sim *n_halfgal * n_day)){
  #Find row in growth parameter data that corresponds to allele sample
  allele_index <- which(growth_import$ST == data3$AT[i]) 
  
  #Calculate the log10N count using our new growth parameters
  #multiply log10N(data$day[i] by 24 since I did it in hours, not days, and ariel did it in days
  data3$count[i] <- log10N(data3$day[i]*24, data3$newLag[i], data3$newMu[i],data3$logMPN_init[i],growth_import$Nmax[allele_index],growth_import$model_name[allele_index])

  #log info
  if(i %% 1000==0) {
    cat(strftime(Sys.time(),"%Y-%m-%d %H:%M:%S",'\n'), ' Iteration ', i, "complete",'/n')
  }
}


save.image(file = "TempNeg60_ppcmodel_1.11init.RData")

