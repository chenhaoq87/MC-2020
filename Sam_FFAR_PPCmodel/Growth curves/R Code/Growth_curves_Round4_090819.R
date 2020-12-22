#growth curves round 4
#edited by SL on 9/08/19

library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
library(AICcmodavg)
#insert growth model data
growthmodel2 <- read.csv("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Growth curves/R Code/Growth_Curve_Round4_RAWDATA_SL_090819.csv")
names(growthmodel2) <- c("t","Isolate","LOG10N","Count")

R101587 <- growthmodel2[growthmodel2$Isolate=="FSL R10-1587", ]

#Run Buchanan growth model
R101587modelb <- nls(buchanan, R101587,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))
                   # control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))

#Run Gompertz growth model
R101587modelg <- nls(gompertzm, R101587,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R101587modelba <- nls(baranyi, R101587,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#output Bic table for each model
candidate_models <- list()
candidate_models[[1]] <- R101587modelb
candidate_models[[2]] <- R101587modelg
candidate_models[[3]] <- R101587modelba
mod.names <- c("Buchanan", "Gompertz", "Baryani")
title_string <- paste("Isolate R10-1587")
title_string
output_bic <- bictab(cand.set = candidate_models, modnames = mod.names, sort = TRUE)
print(title_string)
print(output_bic)

#the next two lines say to save the plot as a png file instead of displaying it
new_filename = paste("Plot", "isolate R10-1587",sep="_")
png(file = paste(new_filename, ".png", sep=''), width=1600, height = 900)
#draw the plot
plot(R100056$t, R100056$LOG10N, xlim = c(0,240), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100056$t),max(R100056$t)))
lines(new$t,predict(R101587modelg,newdata=new), col="red")
lines(new$t,predict(R101587modelb,newdata = new), col="deepskyblue")
lines(new$t,predict(R101587modelba,newdata = new))
#add a nice title containing the isolate name and bioreplic
#the paste function lets you add data into a string
title_string <- paste("Isolate R10-1587")
s <- paste0("Buchanan = Blue, Baranyi = Black, Gompertz = Red")
title(main = title_string, sub = s)
dev.off()

R101587buch<-coef(R101587modelb)
R101587gomp<-coef(R101587modelg)
R101587bary<-coef(R101587modelba)


R100587 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0587", ]

#Run Buchanan growth model
#R100084modelb <- nls(buchanan, R100084, 
                     #list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                    # control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))
#Run Gompertz growth model
R100587modelb <- nls(buchanan, R100587,
                     list(lag = 4, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Gompertz growth model
R100587modelg <- nls(gompertzm, R100587,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R100587modelba <- nls(baranyi, R100587,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))


#output Bic table for each model
candidate_modelsb <- list()
candidate_modelsb[[1]] <- R100587modelb
candidate_modelsb[[2]] <- R100587modelg
candidate_modelsb[[3]] <- R100587modelba
mod.namesb <- c("Buchanan", "Gompertz", "Baryani")
title_stringb <- paste("Isolate R10-0587")
title_stringb
output_bicb <- bictab(cand.set = candidate_modelsb, modnames = mod.names, sort = TRUE)
print(title_stringb)
print(output_bicb)


#the next two lines say to save the plot as a png file instead of displaying it
new_filename = paste("Plot", "isolate R10-0587",sep="_")
png(file = paste(new_filename, ".png", sep=''), width=1600, height = 900)
#draw the plot
plot(R100587$t, R100587$LOG10N, xlim = c(0,240), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100587$t),max(R100587$t)))
lines(new$t,predict(R100587modelg,newdata=new), col="red")
lines(new$t,predict(R100587modelb,newdata = new), col="deepskyblue")
lines(new$t,predict(R100587modelba,newdata = new))
#add a nice title containing the isolate name and bioreplic
#the paste function lets you add data into a string
title_string <- paste("Isolate R10-0587")
s <- paste0("Buchanan = Blue, Baranyi = Black, Gompertz = Red")
title(main = title_string, sub = s)
dev.off()


R100587buch<-coef(R100587modelb)
R100587gomp<-coef(R100587modelg)
R100587bary<-coef(R100587modelba)




Growthcurvestroubleshooting<- rbind(R101587buch, R101587gomp, R101587bary,R100587buch, R100587gomp, R100587bary)
write.csv(Growthcurvestroubleshooting, "Growth_curves_round_four_090819.csv")
