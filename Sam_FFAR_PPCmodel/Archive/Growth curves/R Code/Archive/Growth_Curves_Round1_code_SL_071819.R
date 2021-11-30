library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
#insert growth model data
growthmodel2 <- read.csv("Growth_Curve_Round1_SL_071219.csv")
names(growthmodel2) <- c("t","Isolate","LOG10N","Count", "BioRep")

R100084 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0084", ]

#Run Buchanan growth model
R100084modelb <- nls(buchanan, R100084,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100084modelg <- nls(gompertzm, R100084,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R100084modelba <- nls(baranyi, R100084,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

BIC(R100084modelb)
BIC(R100084modelg)
BIC(R100084modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100084$t, R100084$LOG10N, xlim = c(0,360), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100084$t),max(R100084$t)))
lines(new$t,predict(R100084modelb,newdata=new))

overview(R100084modelb)
R100084bioone<-coef(R100084modelb)

R100990 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0990", ]

#Run Buchanan growth model
R100990modelb <- nls(buchanan, R100990,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100990modelg <- nls(gompertzm, R100990,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R100990modelba <- nls(baranyi, R100990,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

BIC(R100990modelb)
BIC(R100990modelg)
BIC(R100990modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100990$t, R100990$LOG10N, xlim = c(0,360), ylim = c(3,8), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100990$t),max(R100990$t)))
lines(new$t,predict(R100990modelb,newdata=new))

overview(R100990modelb)
R100990bioone<-coef(R100990modelb)

#parametersR100990<-overview(R100990modelg)
#write.csv(parametersR100990, "parametersR100990_071819.csv")

R103286 <- growthmodel2[growthmodel2$Isolate=="FSL R10-3286", ]

#Run Buchanan growth model
R103286modelb <- nls(buchanan, R103286,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R103286modelg <- nls(gompertzm, R103286,
                    list(lag =6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R103286modelba <- nls(baranyi, R103286,
                     list(lag = 2, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

BIC(R103286modelb)
BIC(R103286modelg)
BIC(R103286modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R103286$t, R103286$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R103286$t),max(R103286$t)))
lines(new$t,predict(R103286modelb,newdata=new))

overview(R103286modelb)
R103286bioone<-coef(R103286modelb)

R102381 <- growthmodel2[growthmodel2$Isolate=="FSL R10-2381", ]

#Run Buchanan growth model
R102381modelb <- nls(buchanan, R102381,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R102381modelg <- nls(gompertzm, R102381,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R102381modelba <- nls(baranyi, R102381,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

BIC(R102381modelb)
BIC(R102381modelg)
BIC(R102381modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Buchanan has lower BIC value
#a rough plot to see how our model looks
plot(R102381$t, R102381$LOG10N, xlim = c(0,360), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R102381$t),max(R102381$t)))
lines(new$t,predict(R102381modelb,newdata=new))

overview(R102381modelb)
R102381bioone<-coef(R102381modelb)

rbind(R100084bioone, R100990bioone, R103286bioone, R102381bioone)

