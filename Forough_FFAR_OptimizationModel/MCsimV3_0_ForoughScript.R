##################################################
#  Monte Carlo Psychrotolerant Sporeformer Simulation v3.0
##################################################
## Project: Dairy Spoilage Model
## Script purpose: Simulate spore growth for half gallon samples of milk over a given
##                 number of days.
## Date:  October 01, 2019
## Original Author: Ariel Buehler, Cornell University, ajb466@cornell.edu
##
## V3.0 by: Mike Phillips, Cornell University, mdp38@cornell.edu
##################################################
## Notes: Version  3.0 adds a series of stages to the temperature
## This can be used to model things like variation between or among suppy chains. 
##################################################
# Purpose: This script uses the mc sim v3.0 to calculate the percentage spoiled per day
#          given an initial count.
# How to use: Modify line 202 with your init counts (given in log values) and then run.
#             A plot of percent spoiled per day is given and the first day spoiled is
#             printed to the screen. 
#             To change what is calculated or displayed look from 205 to the
#             end of the file, you should not need to change any other structural things.
library(ggplot2)
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path ))
print( getwd() )


#use a seed value for reproducibility
seed_value = 42;
set.seed(seed_value)


# Utility Functions   -----
##muAtNewTemp
#Purpose: Calculate the new mu parameter at new temperature.
#Params:  newTemp: the new temperature for which we calculate mu
#         oldMu: the previous mu value to adjust
#         oldTemp: the temperature corresponding to previous mu
#         T0:    Parameter used to calculate new mu
muAtNewTemp <- function(newTemp, oldMu, oldTemp = 6, T0 = -3.62) {
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  
  return(newMu)
}

##lagAtNewTemp
#Purpose: Calculate the new lag parameter at new temperature.
#Params:  newTemp: the new temperature for which we calculate lag
#         oldLag: the previous lag value to adjust
#         oldTemp: the temperature corresponding to previous lag
#         T0:    Parameter used to calculate new lag
lagAtNewTemp <- function (newTemp, oldLag, oldTemp = 6, T0 = -3.62) {
  numerator <- oldTemp -T0
  denom <- newTemp - T0
  newLag <- ( (numerator / denom)^2) * oldLag
  return(newLag)
}
#Growth Models
buchanan_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax){
  ans <- LOG10N0 + (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) *     log(10)/mumax)) * mumax * (t - lag)/log(10) + (t >= lag) * (t > (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * (LOG10Nmax -     LOG10N0)
  return(ans)
}
gompertz_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10N0 + (LOG10Nmax - LOG10N0) * exp(-exp(mumax * exp(1) * 
                                                      (lag - t)/((LOG10Nmax - LOG10N0) * log(10)) + 1))
  return(ans)
}
baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * 
                                                          t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax - 
                                                                                                             LOG10N0)))
  return(ans)
}

#Function to calculate log10N
#Wwrapper function because it calls the proper model
#Purpose: This implements the growth model
log10N_func <- function(t, lag, mumax, LOG10N0, LOG10Nmax, model_name="buchanan") {
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



frequency_file <- "Frequency.csv"
growth_file <- "GrowthParameters.csv"
n_halfgal <-10
n_day <- 24
start_day <- 1

#runSim
#Purpose: Run the simulation with the given values
runSim <- function( n_sim = 1000, initMPN = 1.70, 
                    seed_value = 42) {
  # Data frame creation and setup   ----
  #Set up data frame to store count at each day
  #Size is for n_sim bulk tanks, n_half_gal half gallon lots, n_day days

  
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
  


  params <- list ("seed" = seed_value, "initMPN" = initMPN,
                  "freq_file"=frequency_file,
                  "growth_file" = growth_file,
                  "bulk_tanks" = n_sim, "half_gallon_lots"= n_halfgal, "days" =n_day
                  )
  
  
  #Import frequency data and get the rpoB allelic type
  freq_import <- read.csv(frequency_file, stringsAsFactors = FALSE, header = TRUE)
  freq_data = freq_import$rpoB.allelic.type
  
  #Import growth parameter data
  growth_import <-read.csv(growth_file, stringsAsFactors = FALSE)
  

  # Calculate samples used in the monte carlo   ----
  
  #
  #MPN For this version is given by a constant
  data$logMPN_init <- initMPN #Add initial logMPN to data frame

  #Temperature data
  stages <- read.csv("temp_stages.csv", stringsAsFactors = F, comment.char = "#")
  
  #Generate initial MPN for each half gallon from Poisson distribution
  #Also sample AT for each half gallon
  allele <- vector()
  temps <- vector()
  for (i in 1:n_sim){
    allele_samp <- rep(sample(freq_data, n_halfgal, replace = T), times = n_day)
    allele <- c(allele, allele_samp)
    #now calculate temp
    for (j in 1:nrow(stages)){
      stage_row <- stages[j, ]
      n_times <- stage_row$endTime - stage_row$beginTime + 1
      params <- as.numeric(unlist(strsplit(stage_row$parameters, " ")))
      temp_mean <- params[[1]]
      temp_sd <- params[[2]]
      temp_sample <- rep(rnorm(n_halfgal, temp_mean, temp_sd), times = n_times)
      temps <- c(temps, temp_sample)
    }
  }
  #add in temperature
  data$temp <- temps

  #Now we add in those calculations to our original dataframe
  #MPN For this version is given by a constant
  data$logMPN_init <- initMPN 
  data$AT<-allele #Add in AT data
  
  data$AT<-allele #Add in AT data
  
  ##Now we will calculate the log10N for each row in the data frame
  ##Get the AT and day from the data frame, get growth parameters depending on the AT
  # Simulation   ----
  for (i in 1:(n_sim *n_halfgal * n_day)){
    #Find row in growth parameter data that corresponds to allele sample
    allele_index <- which(growth_import$rpoBAT == data$AT[i]) 
    
    #calculate the new growth parameters using the square root model and our
    #sampled temperature
    newT <- data$temp[i]
    newLag <- lagAtNewTemp(newT, growth_import$lag[allele_index])
    newMu <-  muAtNewTemp(newT, growth_import$mumax[allele_index])
    
    #Calculate the log10N count using our new growth parameters
    data$count[i] <- log10N_func(data$day[i], newLag, newMu, data$logMPN_init[i],growth_import$LOG10Nmax[allele_index])
    
  }
  
  retValue = list("results" = data, "params" = params)
  return(retValue)
}


values = c(1.70, 2.24, 2.56, 2.87, 3.30)
print ("Starting sim...")

for(v in values) {
  result <- runSim(n_sim = 100, initMPN = v, seed_value = 42)[[1]];
  
  #get percent reached 5% spoilage per day (defined as count over 6 log)
  days <- seq(1, 24, 1)
  spoiled <- vector(mode = "numeric", length(days))
  spoiled_by_day <- data.frame(days=as.factor(days), spoiled)
  
  for (d in days) {
    numerator <-   length(which(subset(result, result$day == d)$count>=6))
    denominator <- length(subset(result, result$day == d)$count)
    spoiled_by_day[spoiled_by_day$days==d,]$spoiled <- (numerator / denominator) * 100
  }
  
  first_day <- spoiled_by_day[spoiled_by_day$spoiled>5, 'days'][1]
  print(paste0("Init Value: ", v, " reaches 5% spoiled on day ", first_day  ))
  
  ##plot the results
  fig_SBD <- ggplot() + geom_bar(aes(y = spoiled_by_day$spoiled, x = spoiled_by_day$days), data = spoiled_by_day, 
                                 stat = "identity", color = "black", position = "dodge") +
    geom_text(data = spoiled_by_day, 
              aes(x = spoiled_by_day$days, y = spoiled_by_day$spoiled, 
                  label = round(spoiled_by_day$spoiled, 2)), 
              position=position_dodge(width = 1), 
              vjust=-0.25, size = 4) +
    theme_minimal() + labs(y = "Percentage > 6log", x = "Day",
                           title = "Percent of samples spoiled by day",
                           subtitle = paste0("Initial concentration = ", v, " log."))+
    scale_x_discrete(breaks = seq(1, 24, 1), labels=seq(1, 24, 1)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0))
  
  print(fig_SBD)
}




