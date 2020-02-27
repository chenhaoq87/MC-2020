#Consumer Exposure to Mold in Yogurt Predictive Model
#Written by Ariel Buehler, Cornell University, ajb466@cornell.edu

setwd("/Users/arielbuehler/Documents/Documents/Cornell/Research/")
library(triangle)
library(dplyr)
library(truncnorm) #package for truncated normal distribution
library(fitdistrplus)

set.seed(105)
nsim=100
run_sim_domestic <- function(nsim) {
  data.frame(contam=rbinom(n=nsim,size=1,prob=(1/1800)),
             Time_control=rtriangle(n=nsim, a=192, b=792, c=240) + rtriangle(n=nsim, a=2, b=240, c=24),
             T_domestic=rtruncnorm(n=nsim, a=0, b=(17.22), mean=(3.262108), sd=2.623831),
             Time_domestic=rtriangle(n=nsim, a=24, b=240, c=120),
             muopt=rnorm(n=nsim, mean=0.261, sd=0), #sd=0.004
             T_control=runif(n=nsim, min=1, max=5),
             Tmax_mu=rnorm(n=nsim, mean=29.8, sd=0), #sd=0.2
             Tmin_mu=rnorm(n=nsim, mean=-7.6, sd=0), #sd=0.6
             Topt_mu=rnorm(n=nsim, mean= 19.5, sd=0), #sd=0.3
             lagopt=rnorm(n=nsim, mean=0.042, sd=0), #sd=0.002
             Tmax_lag=rnorm(n=nsim, mean =30.1, sd=0), #sd=0.9
             Tmin_lag=rnorm(n=nsim, mean= -6.3, sd=0), #sd=0.8
             Topt_lag=rnorm(n=nsim, mean=23.8, sd=0)) %>% 
      mutate(mu_control=(muopt*(T_control-Tmax_mu)*((T_control-Tmin_mu)^2)) / ((Topt_mu-Tmin_mu)*(((Topt_mu-Tmin_mu)*(T_control-Topt_mu))-((Topt_mu-Tmax_mu)*(Topt_mu+Tmin_mu-2*T_control))))) %>%
      mutate(mu_domestic=(muopt*(T_domestic-Tmax_mu)*((T_domestic-Tmin_mu)^2)) / ((Topt_mu-Tmin_mu)*(((Topt_mu-Tmin_mu)*(T_domestic-Topt_mu))-((Topt_mu-Tmax_mu)*(Topt_mu+Tmin_mu-2*T_domestic))))) %>%
      mutate(lag_control=((Topt_lag-Tmin_lag)*(((Topt_lag-Tmin_lag)*(T_control-Topt_lag))-((Topt_lag-Tmax_lag)*(Topt_lag+Tmin_lag-2*T_control))))/
               (lagopt*(T_control-Tmax_lag)*((T_control-Tmin_lag)^2))) %>%
      mutate(lag_domestic=((Topt_lag-Tmin_lag)*(((Topt_lag-Tmin_lag)*(T_domestic-Topt_lag))-((Topt_lag-Tmax_lag)*(Topt_lag+Tmin_lag-2*T_domestic))))/
               (lagopt*(T_domestic-Tmax_lag)*((T_domestic-Tmin_lag)^2))) %>%
      mutate(Diameter_control=mu_control*(Time_control-lag_control)) %>%
      mutate(Diameter_control_adj=ifelse(Diameter_control<0, 0, Diameter_control)) %>%
      mutate(Diameter_domestic=mu_domestic*(Time_domestic-lag_domestic)) %>%
      mutate(Diameter_domestic_adj=ifelse(Diameter_domestic<0, 0, Diameter_domestic)) %>%
      mutate(Diameter_sum=Diameter_control_adj+Diameter_domestic_adj) %>%
      mutate(Diameter_contam=Diameter_sum*contam) %>%
      mutate(Seemold= 3 < Diameter_contam) -> mold_df 
    
  return(table(mold_df$Seemold)[2]*100/nsim)
}


run100_domestic <- sapply(rep(1e6,100),FUN=run_sim_domestic)
hist(run100_domestic, 
     main="Consumer Exposures to Mold", xlab="Percent of cups with visible mold", 
     col="grey", breaks=10, ylim=c(0, 50), xlim=c(0.045, 0.065))


#What-if the yogurt was regionally distributed??
set.seed(105)
nsim=100
run_sim_small <- function(nsim) {
  data.frame(contam=rbinom(n=nsim,size=1,prob=(1/1800)),
             Time_control=rtriangle(n=nsim, a=24, b=168, c=120) + rtriangle(n=nsim, a=24, b=672, c=336),
             T_domestic=rtruncnorm(n=nsim, a=0, b=(17.22), mean=(3.262108), sd=2.623831),
             Time_domestic=rtriangle(n=nsim, a=24, b=240, c=120),
             muopt=rnorm(n=nsim, mean=0.261, sd=0), #sd=0.004
             T_control=runif(n=nsim, min=1, max=5),
             Tmax_mu=rnorm(n=nsim, mean=29.8, sd=0), #sd=0.2
             Tmin_mu=rnorm(n=nsim, mean=-7.6, sd=0), #sd=0.6
             Topt_mu=rnorm(n=nsim, mean= 19.5, sd=0), #sd=0.3
             lagopt=rnorm(n=nsim, mean=0.042, sd=0), #sd=0.002
             Tmax_lag=rnorm(n=nsim, mean =30.1, sd=0), #sd=0.9
             Tmin_lag=rnorm(n=nsim, mean= -6.3, sd=0), #sd=0.8
             Topt_lag=rnorm(n=nsim, mean=23.8, sd=0)) %>% 
    mutate(mu_control=(muopt*(T_control-Tmax_mu)*((T_control-Tmin_mu)^2)) / ((Topt_mu-Tmin_mu)*(((Topt_mu-Tmin_mu)*(T_control-Topt_mu))-((Topt_mu-Tmax_mu)*(Topt_mu+Tmin_mu-2*T_control))))) %>%
    mutate(mu_domestic=(muopt*(T_domestic-Tmax_mu)*((T_domestic-Tmin_mu)^2)) / ((Topt_mu-Tmin_mu)*(((Topt_mu-Tmin_mu)*(T_domestic-Topt_mu))-((Topt_mu-Tmax_mu)*(Topt_mu+Tmin_mu-2*T_domestic))))) %>%
    mutate(lag_control=((Topt_lag-Tmin_lag)*(((Topt_lag-Tmin_lag)*(T_control-Topt_lag))-((Topt_lag-Tmax_lag)*(Topt_lag+Tmin_lag-2*T_control))))/
             (lagopt*(T_control-Tmax_lag)*((T_control-Tmin_lag)^2))) %>%
    mutate(lag_domestic=((Topt_lag-Tmin_lag)*(((Topt_lag-Tmin_lag)*(T_domestic-Topt_lag))-((Topt_lag-Tmax_lag)*(Topt_lag+Tmin_lag-2*T_domestic))))/
             (lagopt*(T_domestic-Tmax_lag)*((T_domestic-Tmin_lag)^2))) %>%
    mutate(Diameter_control=mu_control*(Time_control-lag_control)) %>%
    mutate(Diameter_control_adj=ifelse(Diameter_control<0, 0, Diameter_control)) %>%
    mutate(Diameter_domestic=mu_domestic*(Time_domestic-lag_domestic)) %>%
    mutate(Diameter_domestic_adj=ifelse(Diameter_domestic<0, 0, Diameter_domestic)) %>%
    mutate(Diameter_sum=Diameter_control_adj+Diameter_domestic_adj) %>%
    mutate(Diameter_contam=Diameter_sum*contam) %>%
    mutate(Seemold= 3 < Diameter_contam) -> mold_df 
  
  return(table(mold_df$Seemold)[2]*100/nsim)
}
run100_small <- sapply(rep(1e6,100),FUN=run_sim_small)

##What-if the domestic fridge storage temp distribution was truncated to 4C?
set.seed(105)
nsim=100
run_sim_fridge <- function(nsim) {
  data.frame(contam=rbinom(n=nsim,size=1,prob=(1/1800)),
             Time_control=rtriangle(n=nsim, a=192, b=792, c=240) + rtriangle(n=nsim, a=2, b=240, c=24),
             T_domestic=rtruncnorm(n=nsim, a=0, b=(17.22), mean=(3.262108-2), sd=2.623831),
             Time_domestic=rtriangle(n=nsim, a=24, b=240, c=120),
             muopt=rnorm(n=nsim, mean=0.261, sd=0), #sd=0.004
             T_control=runif(n=nsim, min=1, max=3),
             Tmax_mu=rnorm(n=nsim, mean=29.8, sd=0), #sd=0.2
             Tmin_mu=rnorm(n=nsim, mean=-7.6, sd=0), #sd=0.6
             Topt_mu=rnorm(n=nsim, mean= 19.5, sd=0), #sd=0.3
             lagopt=rnorm(n=nsim, mean=0.042, sd=0), #sd=0.002
             Tmax_lag=rnorm(n=nsim, mean =30.1, sd=0), #sd=0.9
             Tmin_lag=rnorm(n=nsim, mean= -6.3, sd=0), #sd=0.8
             Topt_lag=rnorm(n=nsim, mean=23.8, sd=0)) %>% 
    mutate(mu_control=(muopt*(T_control-Tmax_mu)*((T_control-Tmin_mu)^2)) / ((Topt_mu-Tmin_mu)*(((Topt_mu-Tmin_mu)*(T_control-Topt_mu))-((Topt_mu-Tmax_mu)*(Topt_mu+Tmin_mu-2*T_control))))) %>%
    mutate(mu_domestic=(muopt*(T_domestic-Tmax_mu)*((T_domestic-Tmin_mu)^2)) / ((Topt_mu-Tmin_mu)*(((Topt_mu-Tmin_mu)*(T_domestic-Topt_mu))-((Topt_mu-Tmax_mu)*(Topt_mu+Tmin_mu-2*T_domestic))))) %>%
    mutate(lag_control=((Topt_lag-Tmin_lag)*(((Topt_lag-Tmin_lag)*(T_control-Topt_lag))-((Topt_lag-Tmax_lag)*(Topt_lag+Tmin_lag-2*T_control))))/
             (lagopt*(T_control-Tmax_lag)*((T_control-Tmin_lag)^2))) %>%
    mutate(lag_domestic=((Topt_lag-Tmin_lag)*(((Topt_lag-Tmin_lag)*(T_domestic-Topt_lag))-((Topt_lag-Tmax_lag)*(Topt_lag+Tmin_lag-2*T_domestic))))/
             (lagopt*(T_domestic-Tmax_lag)*((T_domestic-Tmin_lag)^2))) %>%
    mutate(Diameter_control=mu_control*(Time_control-lag_control)) %>%
    mutate(Diameter_control_adj=ifelse(Diameter_control<0, 0, Diameter_control)) %>%
    mutate(Diameter_domestic=mu_domestic*(Time_domestic-lag_domestic)) %>%
    mutate(Diameter_domestic_adj=ifelse(Diameter_domestic<0, 0, Diameter_domestic)) %>%
    mutate(Diameter_sum=Diameter_control_adj+Diameter_domestic_adj) %>%
    mutate(Diameter_contam=Diameter_sum*contam) %>%
    mutate(Seemold= 3 < Diameter_contam) -> mold_df 
  
  return(table(mold_df$Seemold)[2]*100/nsim)
}


run100_fridge <- sapply(rep(1e6,100),FUN=run_sim_fridge)


#Two histograms overlapping
hist(run100_domestic, main="Consumer Exposures to Mold", xlab="Percent of cups with visible mold", col="darkslategrey", xlim=c(0.040,0.065), ylim=c(0, 50)) 
  hist(run100_fridge, col=rgb(190, 190, 190, max = 255, alpha = 125, names = "gray50"), add=TRUE)
  legend(0.040, 49, c("Lower storage temperature", "Baseline temperature"), pch = c(15, 15), col=c("gray", "darkslategrey"))


  legend(x=0.040, y=32, c("Regionally","baseline"), pch = c(15, 15), col=c("gray", "darkslategrey"))
  
  
  
  ##Histograms with black and white colors
  hist(run100_domestic, main="Consumer Exposures to Mold", xlab="Percent of cups with visible mold", col="black", xlim=c(0.045,0.065), ylim=c(0, 50)) 
  hist(run100_fridge, col=rgb(255, 255, 255, max = 255, alpha = 125, names = "gray50"), add=TRUE)
  legend(0.040, 49, c("Lower storage temperature", "Baseline temperature"), pch = c(15, 15), col=c("whitesmoke", "black"))

  
 
    
#chobani probability of contaminated =10/22,000
#Chobani time from processing plant to retail store: in hours: rltriangle(a=192, b=792, c=240)
#Chobani time from retail to consumer: in hours: rltriangle(n=nsim, a=2, b=240, c=24)
#Chobani time in consumer fridge in hours: rltriangle(n=nsim, a=24, b=240, c=120)
#Ithaca Milk time from processing plant to retail store: Mean=5days, Min=1 day, Max=7days (rltriangle a=1, b=7, c=5) in hours: (rltriangle a=24, b=168, c=120)
#Ithaca Milk time on retail store shelf before being bought by consumer: Mean=2 weeks, min=1 day, max=4weeks (rltriangle a=1, b=28, c=14) in hours: (rltriangle a=24, b=672, c=336)


#Storage temperature in domestic fridge, degrees C, using EcoSure 2007 study
Yogurt_domestic_import <- read.csv("/Users/arielbuehler/Documents/Documents/Cornell/Research/Chobani Mold Model/Yogurt_domestictemps.csv")
Yogurt_domestic_degC <- Yogurt_domestic_import$DegC
temp_dist <- fitdist(Yogurt_domestic_degC, "norm")

curve(dnorm(x, mean=3.26, sd=2.62), add=TRUE)
hist(Yogurt_domestic_degC, xlab="Temperature (C) of yogurt in domestic refrigerator")

sampletemp = function(n) { 
  sample(Yogurt_domestic_degC, n, replace = T) 
}

#Growth (mm) during truck and retail, using Cardinal Model with inflection developed by Rosso et al. 1999
#Penicillium commune values from Gougouli, 2011
#tmin = -7.6 +- 0.6
#tmax = 29.8 +-0.2
#topt = 19.5 +- 0.3
#muopt = 0.257 +- 0.004
mu_func = function(muopt,t,tmax_mu,tmin_mu,topt_mu){
  ans <-(muopt*(t-tmax_mu)*((t-tmin_mu)^2)) / ((topt_mu-tmin_mu)*(((topt_mu-tmin_mu)*(t-topt_mu))-((topt_mu-tmax_mu)*(topt_mu+tmin_mu-2*t))))
  return(ans)
}

#Lag (hrs) during truck and retail, using Cardinal Model with inflection developed by Rosso et al. 1999 and Penicillium commune values from Gougouli, 2011
#tmin=-6.3 +-0.8
#tmax=30.1 +-0.9
#topt = 23.8 +-0.4
#lagopt=0.040 +- 0.002
lag_func = function(lagopt,t,tmax_lag,tmin_lag,topt_lag){
  ans <- ((topt_lag-tmin_lag)*(((topt_lag-tmin_lag)*(t-topt_lag))-((topt_lag-tmax_lag)*(topt_lag+tmin_lag-2*t))))/
  (lagopt*(t-tmax_lag)*((t-tmin_lag)^2))
  return(ans)
}


#this answer is how many cups could have mold by the time a consumer opens cup.
#Feed contaminated number into the time_sum vector to see if the cup reaches consumer before or after mold growth
#If cup reaches consumer after mold growth =1
#If cup reaches consumer before mold growth =0

######################################################################################
  
#nsim=1e6

#run_sim <- function(nsim) {
#  data.frame(tvg=rltriangle(n=nsim, a=9, b=34, c=23),
 #            contam=rbinom(n=nsim,size=1,prob=(10/22000)),
  #           time_sum=rltriangle(n=nsim, a=8, b=33, c=10) +
   #            rltriangle(n=nsim, a=0.083, b=10, c=1) +
    #           rltriangle(n=nsim, a=1, b=10, c=5)) %>%
    #mutate(combined=tvg/contam) %>%
    #mutate(consumer_sees_mold=combined<time_sum) -> sim_df
  
#  return(table(sim_df$consumer_sees_mold)[2]*100/nsim)
#}

#run100 <- sapply(rep(1e6,100),FUN=run_sim)
#hist(run100, main="Consumer Exposures to Yeast and Mold", xlab="Percent of cups with Mold")


#If time_sum is less than time to visible growth, consumer eats yogurt with no yeast/mold contamination.
#If time_sum is more than time to visible growth, consumer opens yogurt WITH yeast/mold contamination


#complaints per lot