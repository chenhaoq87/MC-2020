#growth curves round 4
#edited by SL on 9/08/19

install.packages("minpack.lm")
library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
library(AICcmodavg)
library(minpack.lm)
#insert growth model data
#growthmodel2 <- read.csv("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Growth curves/R Code/Growth_Curve_Round4_RAWDATA_SL_090819.csv")
growthmodel2 <- read.csv("Growth_Curve_Round4_RAWDATA_SL_090819.csv")
names(growthmodel2) <- c("t","Isolate","LOG10N","Count", "AVG")

R101587 <- growthmodel2[growthmodel2$Isolate=="FSL R10-1587", ]

#Run Buchanan growth model
R101587modelb <- nlsLM(start=c(lag = 0, mumax=0.1, LOG10N0=3, LOG10Nmax=8),
                        formula = buchanan, 
                        data = R101587,
                        lower=c(lag = 0, mumax=0, LOG10N0=0, LOG10Nmax=0))
                   # control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))

#Run Gompertz growth model
R101587modelg <- nlsLM(formula = gompertzm, 
                        data = R101587,
                        start = c(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                        lower=c(lag = 0, mumax=0, LOG10N0=0, LOG10Nmax=0))

#Run Baranyi growth model
R101587modelba <- nlsLM(formula = baranyi, 
                      data = R101587,
                     start = c(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                     lower=c(lag = 0, mumax=0, LOG10N0=0, LOG10Nmax=0))

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
plot(R101587$t, R101587$LOG10N, xlim = c(0,240), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R101587$t),max(R101587$t)))
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
#CANNOT FIT BUCHANAN TO THIS DATA
#R100587modelb <- nlsLM(start=c(lag = 199, mumax=0.05, LOG10N0=3, LOG10Nmax=9),
                       #formula = buchanan, 
                       #data = R100587,
                       #lower=c(lag = 0, mumax=0, LOG10N0=0, LOG10Nmax=0))
# control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))

#Run Gompertz growth model
R100587modelg <- nlsLM(formula = gompertzm, 
                       data = R100587,
                       start = c(lag = 100, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                       lower=c(lag = 0, mumax=0, LOG10N0=0, LOG10Nmax=0))

#Run Baranyi growth model
R100587modelba <- nlsLM(formula = baranyi, 
                        data = R100587,
                        start = c(lag = 125, mumax=0.1, LOG10N0=1, LOG10Nmax=8), 
                        lower=c(lag = 0, mumax=0, LOG10N0=0, LOG10Nmax=0))


#output Bic table for each model
#REMOVE BUCHANAN AS IT'S NOT VALID FOR THIS SET
candidate_modelsb <- list()
#candidate_modelsb[[1]] <- R100587modelb
candidate_modelsb[[1]] <- R100587modelg
candidate_modelsb[[2]] <- R100587modelba
mod.namesb <- c("Gompertz", "Baryani")
title_stringb <- paste("Isolate R10-0587")
title_stringb
output_bicb <- bictab(cand.set = candidate_modelsb, modnames = mod.namesb, sort = TRUE)
print(title_stringb)
print(output_bicb)


#the next two lines say to save the plot as a png file instead of displaying it
new_filename = paste("Plot", "isolate R10-0587",sep="_")
png(file = paste(new_filename, ".png", sep=''), width=1600, height = 900)
#draw the plot
plot(R100587$t, R100587$LOG10N, xlim = c(0,320), ylim = c(3,12), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100587$t),max(R100587$t)))
lines(new$t,predict(R100587modelg,newdata=new), col="red")
#lines(new$t,predict(R100587modelb,newdata = new), col="deepskyblue")
lines(new$t,predict(R100587modelba,newdata = new))
#add a nice title containing the isolate name and bioreplic
#the paste function lets you add data into a string
title_string <- paste("Isolate R10-0587")
s <- paste0("Baranyi = Black, Gompertz = Red")
title(main = title_string, sub = s)
dev.off()


#R100587buch<-coef(R100587modelb)
R100587gomp<-coef(R100587modelg)
R100587bary<-coef(R100587modelba)




Growthcurvestroubleshooting<- rbind(R101587buch, R101587gomp, R101587bary,R100587gomp, R100587bary)
write.csv(Growthcurvestroubleshooting, "Growth_curves_round_four_091019.csv")
