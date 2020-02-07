## Utility Functions
# Last edited: 020620

## (1) muAtNewTemp

# Source: this function uses Ratkowsky's square root model which describes the effect of temperature on the growth of microorganisms (https://www.ncbi.nlm.nih.gov/pubmed/22417595)

# Purpose: Calculate the new mu parameter at new temperature.

# Parameters: 
# (i) newTemp: New temperature for which we calculate mu
# (ii) oldMu: Previous mu value to adjust
# (iii) oldTemp: Temperature corresponding to previous mu; Here, oldTemp = 6C since growth curve experiments were originally performed at 6C
# (iv) T0: Parameter used to calculate new mu; Note: Here, T0 = -3.62C was determined using Ratkowsky's square root model 
#and Paenibacillus ordorifer growth curves obtained at 4, 7, and 32C in BHI broth (N.H. Martin unpublished data)

# Function:
muAtNewTemp <- function(newTemp, oldMu, oldTemp = 6, T0 = -3.62) { 
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  
  return(newMu)
}


## (2) adjustLag

# Purpose: Adjust the lag phase based on the Zwietering 1994 paper.

# Parameters:  
# (i) t: Current timestep (in days)
# (ii) oldLag: Lag time (in days) at the previous temperature 
# (iii) newLag: Lag time (in days) at the current temparature 
# (iv) restartExp: If true, then lag phase restarts even if already in exponential growth phase
# (v) adjustmentConstant: Amount to adjust lag, the paper recommends 0.25

# Function:
adjustLag <- function (t, oldLag, newLag, restartExp = T, adjustmentConstant = 0.25) {
  #determine the amount of lag phase completed
  remainingLag <- 1 - (t / oldLag)
  if(restartExp) {
    remainingLag <- ifelse(remainingLag < 0, 0, remainingLag)
  }
  else {
    adjustedLag <- ifelse(remainingLag <=0, oldLag,
                          t + remainingLag * newLag + adjustmentConstant*newLag)
  }
  
  adjustedLag <- t + remainingLag*newLag + adjustmentConstant*newLag
  return(adjustedLag)
}


## (3) lagAtNewTemp

# Purpose: Calculate the new lag parameter at new temperature.

# Parameters:
# (i) t: Current timestep (in days)
# (ii) newTemp: New temperature for which we calculate lag
# (iii) oldLag: Previous lag value to adjust
# (iv) oldTemp: Temperature corresponding to previous lag
# (v) T0: Parameter used to calculate new mu; Note: Here, T0 = -3.62C was determined using Ratkowsky's square root model 
#and Paenibacillus ordorifer growth curves obtained at 4, 7, and 32C in BHI broth (N.H. Martin unpublished data)

# Function:
lagAtNewTemp <- function (t, newTemp, oldLag, oldTemp = 6, T0 = -3.62) {
  numerator <- oldTemp -T0
  denom <- newTemp - T0
  newLag <- ( (numerator / denom)^2) * oldLag
  return(newLag)
}

getPrevRow <- function(df, sim_run, milk_unit, day) {
  old_temp <- df[df$bulk_tank == sim_run & df$milk_unit == milk_unit & df$day==day-1,] 
}


### Growth Models ###
# Source: All equations for the 3 models below are copied from the Nlms package in R (https://rdrr.io/cran/nlsMicrobio/src/R/growthmodels.R)


## (4) buchanan_log10N

# Purpose:

# Parameters:

# Function:
buchanan_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax){
  ans <- LOG10N0 + (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * mumax * (t - lag)/log(10) + (t >= lag) * (t > (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * (LOG10Nmax -     LOG10N0)
  return(ans)
}


## (5) gompertz_log10N

# Purpose:

# Parameters:

# Function:
gompertz_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10N0 + (LOG10Nmax - LOG10N0) * exp(-exp(mumax * exp(1) * (lag - t)/((LOG10Nmax - LOG10N0) * log(10)) + 1))
  return(ans)
}


## (6) baranyi_log10N

# Purpose:

# Parameters:

baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax - LOG10N0)))
  return(ans)
}


## (7) log10N_func

#Purpose: calculate log10N and implements the growth model

#Parameters:

# Function:
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