#growth curves round 1 and 3 troubleshooting
#edited by SL on 8/25/19

library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
library(AICcmodavg)
#insert growth model data
growthmodel2 <- read.csv("Growth_curve_round1and3RAWDATAcombined_SL_080619.csv")
names(growthmodel2) <- c("t","Isolate","LOG10N","Count")

R100056 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0056", ]

#Run Buchanan growth model
R100056modelb <- nls(buchanan, R100056,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100056modelg <- nls(gompertzm, R100056,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R100056modelba <- nls(baranyi, R100056,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#output Bic table for each model
candidate_models <- list()
candidate_models[[1]] <- R100056modelb
candidate_models[[2]] <- R100056modelg
candidate_models[[3]] <- R100056modelba
mod.names <- c("Buchanan", "Gompertz", "Baryani")
title_string <- paste("Isolate R10-0056")
title_string
output_bic <- bictab(cand.set = candidate_models, modnames = mod.names, sort = TRUE)
print(title_string)
print(output_bic)

#the next two lines say to save the plot as a png file instead of displaying it
new_filename = paste("Plot", "isolate R10-0056",sep="_")
png(file = paste(new_filename, ".png", sep=''), width=1600, height = 900)
#draw the plot
plot(R100056$t, R100056$LOG10N, xlim = c(0,240), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100056$t),max(R100056$t)))
lines(new$t,predict(R100056modelg,newdata=new), col="red")
lines(new$t,predict(R100056modelb,newdata = new), col="deepskyblue")
lines(new$t,predict(R100056modelba,newdata = new))
#add a nice title containing the isolate name and bioreplic
#the paste function lets you add data into a string
title_string <- paste("Isolate R10-0056")
s <- paste0("Buchanan = Blue, Baranyi = Black, Gompertz = Red")
title(main = title_string, sub = s)
dev.off

R100056buch<-coef(R100056modelb)
R100056gomp<-coef(R100056modelg)
R100056bary<-coef(R100056modelba)


R100084 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0084", ]

#Run Buchanan growth model
#R100084modelb <- nls(buchanan, R100084, 
                     #list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                    # control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))
#Run Gompertz growth model
R100084modelb <- nls(buchanan, R100084,
                     list(lag = 4, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Gompertz growth model
R100084modelg <- nls(gompertzm, R100084,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R100084modelba <- nls(baranyi, R100084,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))


#output Bic table for each model
candidate_modelsb <- list()
candidate_modelsb[[1]] <- R100084modelb
candidate_modelsb[[2]] <- R100084modelg
candidate_modelsb[[3]] <- R100084modelba
mod.namesb <- c("Buchanan", "Gompertz", "Baryani")
title_stringb <- paste("Isolate R10-0084")
title_stringb
output_bicb <- bictab(cand.set = candidate_modelsb, modnames = mod.names, sort = TRUE)
print(title_stringb)
print(output_bicb)


#the next two lines say to save the plot as a png file instead of displaying it
new_filename = paste("Plot", "isolate R10-0084",sep="_")
png(file = paste(new_filename, ".png", sep=''), width=1600, height = 900)
#draw the plot
plot(R100084$t, R100084$LOG10N, xlim = c(0,240), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100084$t),max(R100084$t)))
lines(new$t,predict(R100084modelg,newdata=new), col="red")
lines(new$t,predict(R100084modelb,newdata = new), col="deepskyblue")
lines(new$t,predict(R100084modelba,newdata = new))
#add a nice title containing the isolate name and bioreplic
#the paste function lets you add data into a string
title_string <- paste("Isolate R10-0084")
s <- paste0("Buchanan = Blue, Baranyi = Black, Gompertz = Red")
title(main = title_string, sub = s)
dev.off


R100084buch<-coef(R100084modelb)
R100084gomp<-coef(R100084modelg)
R100084bary<-coef(R100084modelba)


R100908 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0908", ]

#Run Buchanan growth model
R100908modelb <- nls(buchanan, R100908,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100908modelg <- nls(gompertzm, R100908,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R100908modelba <- nls(baranyi, R100908,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#output Bic table for each model
candidate_models <- list()
candidate_models[[1]] <- R100908modelb
candidate_models[[2]] <- R100908modelg
candidate_models[[3]] <- R100908modelba
mod.names <- c("Buchanan", "Gompertz", "Baryani")
title_string <- paste("Isolate R10-0908")
title_string
output_bic <- bictab(cand.set = candidate_models, modnames = mod.names, sort = TRUE)
print(title_string)
print(output_bic)


#the next two lines say to save the plot as a png file instead of displaying it
new_filename = paste("Plot", "isolate R10-0908",sep="_")
png(file = paste(new_filename, ".png", sep=''), width=1600, height = 900)
#draw the plot
plot(R100908$t, R100908$LOG10N, xlim = c(0,240), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100908$t),max(R100908$t)))
lines(new$t,predict(R100908modelg,newdata=new), col="red")
lines(new$t,predict(R100908modelb,newdata = new), col="deepskyblue")
lines(new$t,predict(R100908modelba,newdata = new))
#add a nice title containing the isolate name and bioreplic
#the paste function lets you add data into a string
title_string <- paste("Isolate R10-0908")
s <- paste0("Buchanan = Blue, Baranyi = Black, Gompertz = Red")
title(main = title_string, sub = s)
dev.off


R100908buch<-coef(R100908modelb)
R100908gomp<-coef(R100908modelg)
R100908bary<-coef(R100908modelba)



R101432 <- growthmodel2[growthmodel2$Isolate=="FSL R10-1432", ]

#Run Buchanan growth model
R101432modelb <- nls(buchanan, R101432,
                     list(lag = 3, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Gompertz growth model
R101432modelg <- nls(gompertzm, R101432,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R101432modelba <- nls(baranyi, R101432,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))


#output Bic table for each model
candidate_models <- list()
candidate_models[[1]] <- R101432modelb
candidate_models[[2]] <- R101432modelg
candidate_models[[3]] <- R101432modelba
mod.names <- c("Buchanan", "Gompertz", "Baryani")
title_string <- paste("Isolate R10-1432")
title_string
output_bic <- bictab(cand.set = candidate_models, modnames = mod.names, sort = TRUE)
print(title_string)
print(output_bic)


#the next two lines say to save the plot as a png file instead of displaying it
new_filename = paste("Plot", "isolate R10-1432",sep="_")
png(file = paste(new_filename, ".png", sep=''), width=1600, height = 900)
#draw the plot
plot(R101432$t, R101432$LOG10N, xlim = c(0,240), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R101432$t),max(R101432$t)))
lines(new$t,predict(R101432modelg,newdata=new), col="red")
lines(new$t,predict(R101432modelb,newdata = new), col="deepskyblue")
lines(new$t,predict(R101432modelba,newdata = new))
#add a nice title containing the isolate name and bioreplic
#the paste function lets you add data into a string
title_string <- paste("Isolate R10-1432")
s <- paste0("Buchanan = Blue, Baranyi = Black, Gompertz = Red")
title(main = title_string, sub = s)
dev.off


R101432buch<-coef(R101432modelb)
R101432gomp<-coef(R101432modelg)
R101432bary<-coef(R101432modelba)


#MUST FIX THIS, ERROR MESSAGE
R100701 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0701", ]

#Run Buchanan growth model
R100701modelb <- nls(buchanan, R100701,
                     list(lag = 3, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Gompertz growth model
R100701modelg <- nls(gompertzm, R100701,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R100701modelba <- nls(baranyi, R100701,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))


#output Bic table for each model
candidate_models <- list()
candidate_models[[1]] <- R100701modelb
candidate_models[[2]] <- R100701modelg
candidate_models[[3]] <- R100701modelba
mod.names <- c("Buchanan", "Gompertz", "Baryani")
title_string <- paste("Isolate R10-0701")
title_string
output_bic <- bictab(cand.set = candidate_models, modnames = mod.names, sort = TRUE)
print(title_string)
print(output_bic)


#the next two lines say to save the plot as a png file instead of displaying it
new_filename = paste("Plot", "isolate R10-0701",sep="_")
png(file = paste(new_filename, ".png", sep=''), width=1600, height = 900)
#draw the plot
plot(R100701$t, R100701$LOG10N, xlim = c(0,240), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100701$t),max(R100701$t)))
lines(new$t,predict(R100701modelg,newdata=new), col="red")
lines(new$t,predict(R100701modelb,newdata = new), col="deepskyblue")
lines(new$t,predict(R100701modelba,newdata = new))
#add a nice title containing the isolate name and bioreplic
#the paste function lets you add data into a string
title_string <- paste("Isolate R10-0701")
s <- paste0("Buchanan = Blue, Baranyi = Black, Gompertz = Red")
title(main = title_string, sub = s)
dev.off



R100701buch<-coef(R100701modelb)
R100701gomp<-coef(R100701modelg)
R100701bary<-coef(R100701modelba)


Growthcurvestroubleshooting<- rbind(R100056buch, R100056gomp, R100056bary,R100084buch, R100084gomp, R100084bary,R100908buch, R100908gomp, R100908bary,R101432buch, R101432gomp, R101432bary)
write.csv(Growthcurvestroubleshooting, "Growth_curves_round_oneandthree_troubleshooting_082619.csv")
