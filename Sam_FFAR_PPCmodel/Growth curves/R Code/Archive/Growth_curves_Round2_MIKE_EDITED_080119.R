#Growth curves loop example for Sam
#written by MDP on 7/26/2019

library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
library(AICcmodavg)


#get the data
growthmodel2 <- read.csv("Growth_Curve_Round2_allbio_SL_071919.csv", stringsAsFactors = FALSE)
names(growthmodel2) <- c("t","Isolate","LOG10N","Count", "BioRep") 

#get the list of all the possible isolates
isolate_list <- unique(growthmodel2$Isolate)
# The above works by searching the Isolate column and returning all unique values
# We can now use that as our loop.

#we will create a list to put our results in, after all the loops we will bind them
#and write it to a file
output_list <- list()
#we have to keep track of our list position so we can add the new row to the end
list_index <- 1 ##first element will go in 1, then we will increase so next goes in 2, etc

#loop through each isolate
for (iso in isolate_list) {
  
  #loop through each bio replicate
  for (i in seq(1, 3)){
    
    growth_data <- growthmodel2[growthmodel2$Isolate==iso & growthmodel2$BioRep==i, ]
    
    #Run Buchanan growth model
    growth_buch <- nls(buchanan, growth_data,
                       list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                       control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))
    
    growth_gomp <- nls(gompertzm, growth_data,
                         list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8),
                         control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))
    
    growth_bar <- nls(baranyi, growth_data ,
                          list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                          control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))
    
    #it would not save much effort but we could put the next three lines inside
    # a loop as well. Anytime repeating or copying pasting the same stuff a loop
    # is a possible solution
    ##BIC(growth_buch)
    #BIC(growth_gomp )
    #BIC(growth_bar )
    
    #output aic table for each model
    candidate_models <- list()
    candidate_models[[1]] <- growth_buch
    candidate_models[[2]] <- growth_gomp
    candidate_models[[3]] <- growth_bar
    mod.names <- c("Buchanan", "Gompertz", "Baryani")
    title_string <- paste("Isolate ", iso, "and Bioreplicate #", i, sep=" ")
    title_string
    output_aic <- aictab(cand.set = candidate_models, modnames = mod.names, sort = TRUE)
    print(title_string)
    print(output_aic)
    
    
    #the next two lines say to save the plot as a png file instead of displaying it
    new_filename = paste("Plot", "isolate", iso, "biorep", i,  sep="_")
    png(file = paste(new_filename, ".png", sep=''), width=1600, height = 900)
    #draw the plot
    plot(growth_data$t, growth_data$LOG10N, xlim = c(0,360), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
    new = data.frame(t = seq(min(growth_data$t),max(growth_data$t)))
    lines(new$t, predict(growth_buch, newdata = new), col="red")
    lines(new$t,predict(growth_gomp,newdata=new))
    lines(new$t,predict(growth_bar,newdata = new), col="deepskyblue")
    #add a nice title containing the isolate name and bioreplic
    #the paste function lets you add data into a string
    title_string <- paste("Isolate ", iso, "and Bioreplicate #", i, sep=" ")
    s <- paste0("Buchanan = Red, Baranyi = Blue, Gompertz = Black")
    title(main = title_string, sub = s)
    #this line turns off the png function and saves it to a file.
    dev.off()


    #save the results of this isolate and biorep
    results <- coef(growth_buch)
    results["Isolate"] <- iso
    results["biorep"] <- i
    #now add these results to the end of our list of results
    output_list[[list_index]] <- results
    #remember we need to increase i by 1 so the next time we add our new results to the end
    list_index <- list_index + 1
    
  } #closes inner isolate loop
} #closes outer bioreplicate loop


#Yay! Finally we are done with all our processing.
#we need to bind all the items of our list together, 
#the do.call function applies rbind to each item of the list
output <- do.call(rbind, output_list)

#now we can rearrange the columns if we want
output = output[, c(5, 6, 1, 2, 3, 4)]
# Now we just to write our list to a file
write.table(output, file = "Growth_Params.csv", sep=",", row.names=FALSE)

