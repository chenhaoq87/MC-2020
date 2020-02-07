#Psychrotolerant Sporeformer Predictive Model
#Written by Ariel Buehler, Cornell University, ajb466@cornell.edu

setwd("/Users/arielbuehler/Documents/Documents/Cornell/Research/MonteCarloA") #Insert path to where the files are

#Set up data frame to store count at each day
#Size is for n_sim bulk tanks, n_half_gal half gallon lots, n_day days (14-24)
n_sim <-100
n_halfgal <-10
n_day <- 11

BT <- rep(seq(1,n_sim),each=n_halfgal*n_day)
half_gal <- rep(seq(1,n_halfgal),times=n_day*n_sim)
AT <- vector(,n_sim*n_halfgal*n_day)
day <- rep(rep(seq(14,24),each=n_halfgal),times=n_sim)
count <- vector(,n_sim*n_halfgal*n_day)

data <- data.frame(BT,half_gal,AT,day,count)


#Import frequency data
freq_import <- read.csv("ColdgrowATFreq_CrossSect.csv",stringsAsFactors=FALSE, header=TRUE)
freq_data = freq_import[,2]


#Import growth parameter data
growth_import <-read.csv("GrowthParameters.csv",stringsAsFactors=FALSE)
#growth_import$longermumax <- .6*(growth_import$mumax)


#Import initial count logMPN data
initialcount_import <- read.csv("InitialCountsMPN.csv", stringsAsFactors=FALSE)
initialcount_data = initialcount_import[,3]
initialcountlog_data = initialcount_import[,4]

#Function to sample from frequency data (creates n samples)
samplefreq = function(n) { 
  sample(freq_data, n, replace = T) 
}

#rmultinom(1,1,prob=table(freq_data)/length(freq_data))

#Function to calculate log10N
log10N_func = function(t,lag,mumax,LOG10N0,LOG10Nmax){
  ans <- LOG10N0 + (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) *     log(10)/mumax)) * mumax * (t - lag)/log(10) + (t >= lag) *     (t > (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * (LOG10Nmax -     LOG10N0)
  return(ans)
}

#Fit Log Normal Distribution to MPN data
library(fitdistrplus) #load package where your fitdist function is found
lognormaldist <- fitdist(initialcount_data, "lnorm")
summary(lognormaldist)
plot(lognormaldist)
hist(initialcountlog_data, prob=TRUE)
curve(dnorm(x, mean=mean(initialcountlog_data), sd=sd(initialcountlog_data)), add=TRUE)
hist(initialcount_data, prob=TRUE)


#Sample logMPN from normal distribution
mean <- c(-0.7226627)
sd <- c(.9901429)
logMPN_samp <- rnorm(n_sim, mean, sd)#using cross-sectional data without censoring in middle
MPN_samp <- 10^logMPN_samp
MPN_samp_halfgal <- MPN_samp*1900 #MPN per half gallon


#Generate initial MPN for each half gallon from Poisson distribution
#Also generate AT for each half gallon
MPN_init<-vector()
allele <- vector()
for (i in 1:n_sim){
  MPN_init_samp <-rep(rpois(n_halfgal,MPN_samp_halfgal[i]),times=n_day)
  MPN_init<-c(MPN_init,MPN_init_samp)
  
  allele_samp <-rep(samplefreq(n_halfgal),times=n_day)
  allele <- c(allele,allele_samp)
}

#Convert MPN_init from half-gallon to mLs
MPN_init_mL <- MPN_init/1900

data$logMPN_init <- log10(MPN_init_mL) #Add initial logMPN to data frame
data$AT<-allele #Add in AT data

#Vector of times to plot (times we're interested in)
#t_vec <-seq(14,24,by=1)

##Here is the Monte Carlo part
##Now we will calculate the log10N for each row in the data frame
##Get the AT and day from the data frame, get growth parameters depending on the AT


for (i in 1:(n_sim*n_halfgal*n_day)){
  allele_index <- which(growth_import$rpoBAT == data$AT[i]) #Find row in growth parameter data that corresponds to allele sample
  data$count[i] <- log10N_func(data$day[i],growth_import$lag[allele_index],growth_import$mumax[allele_index],data$logMPN_init[i],growth_import$LOG10Nmax[allele_index]) #Calculate log10N
  #data$count[i] <- log10N_func(data$day[i],growth_import$lag[allele_index],rnorm(mean=growth_import$mumax[allele_index],std=.25),data$logMPN_init[i],growth_import$LOG10Nmax[allele_index]) #Calculate log10N (stochastic parameters)
}

##Plot histogram of data 

histogram <- hist(subset(data,data$day == 24)$count,freq=TRUE, col="grey50", )
histogram$breaks
histogram$count

#hist(log10N[11,],breaks=50, xlab="LogCFU/mL") #histogram at 14 days (change 1 to 1-11 to get time 14-24)


##look at each iteration at a time point
#summary_iteration <- data.frame(logMPN_samp, rpoBAT=allele_samp ,count=log10N[1,])

#which have spoiled?
length(which(subset(data, data$day == 21)$count>4.3))/length(subset(data, data$day == 21)$count)
 #length(which(log10N[4,]<4.3))/length(log10N[4,])

#what is the mean count on D21
mean(subset(data, data$day ==20)$count, na.rm=TRUE)
sd(subset(data, data$day ==24)$count, na.rm=TRUE)
