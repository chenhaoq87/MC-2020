#growth curves round 1 and 3 troubleshooting
#edited by SL on 8/24/19

library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
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

BIC(R100056modelb)
BIC(R100056modelg)
BIC(R100056modelba)

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

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
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

#ERROR IN 0084
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

BIC(R100084modelb)
BIC(R100084modelg)
BIC(R100084modelba)

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

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
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

R100701 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0701", ]

#Run Buchanan growth model
R100701modelb <- nls(buchanan, R100701,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Gompertz growth model
R100701modelg <- nls(gompertzm, R100701,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))
#Run Baranyi growth model
R100701modelba <- nls(baranyi, R100701,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

BIC(R100701modelb)
BIC(R100701modelg)
BIC(R100701modelba)

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

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
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

BIC(R100908modelb)
BIC(R100908modelg)
BIC(R100908modelba)

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

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
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

Growthcurvestroubleshooting<- rbind(R100056b, R100056g, R100056ba,R100084b, R100084g, R100084ba,R100701b, R100701g, R100701ba,)
write.csv(Growthcurvesroundone, "Growth_curves_round_oneandthree_troubleshooting_082419.csv")
