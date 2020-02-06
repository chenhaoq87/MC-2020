## Utility Functions
# Last edited: 020620

## (1) muAtNewTemp

# Source: this function uses the Ratkowsky square model which describes the effect of temperature on the growth of microorganisms (https://www.ncbi.nlm.nih.gov/pubmed/22417595)

# Purpose: Calculate the new mu parameter at new temperature.

# Parameters: 
<<<<<<< HEAD
# (i) newTemp: the new temperature for which we calculate mu
# (ii) oldMu: the previous mu value to adjust
# (iii) oldTemp: the temperature corresponding to previous mu
# (iv) T0: parameter used to calculate new mu; 
=======
# (a) newTemp: the new temperature for which we calculate mu
# (b) oldMu: the previous mu value to adjust
# (c) oldTemp: the temperature corresponding to previous mu
# (d) T0: parameter used to calculate new mu; 
>>>>>>> 95b012a8da7b55bcf1565b3db142aab70f093fe2

# Function:
muAtNewTemp <- function(newTemp, oldMu, oldTemp = 6, T0 = -3.62) { 
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  
  return(newMu)
}

# Note: T0 in the above function is estimated to be -3.62C based on growth curves of
#Paenibacillus ordorifer obtained at 4, 7, and 32C in BHI broth (N.H. Martin unpublished data)


## (2) adjustLag

# Purpose: Adjust the lag phase based on the Zwietering 1994 paper.

# Parameters:  
# (i) t: the current timestep
# (ii) oldLag: the lag time at the previous temperature
# (iii) newLag: the lag time at the current temparature.
# (iv) restartExp: If true then lag phase restarts even if already in exponential growth phase.
# (v) adjustmentConstant: The amount to adjust lag, the paper recommends 0.25

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
# (i) newTemp: the new temperature for which we calculate lag
# (ii) oldLag: the previous lag value to adjust
# (iii) oldTemp: the temperature corresponding to previous lag
# (iv) T0: parameter used to calculate new lag

# Function:
lagAtNewTemp <- function (t, newTemp, oldLag, oldTemp = 6, T0 = -3.62) {
  numerator <- oldTemp -T0
  denom <- newTemp - T0
  newLag <- ( (numerator / denom)^2) * oldLag
  return(newLag)
}

getPrevRow <- function(df, sim_run, half_gallon, day) {
  old_temp <- df[df$BT == sim_run & df$half_gal == half_gallon & df$day==day-1,] 
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