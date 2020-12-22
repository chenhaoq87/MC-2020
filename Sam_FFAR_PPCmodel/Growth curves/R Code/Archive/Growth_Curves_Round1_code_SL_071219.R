library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
#insert growth model data
growthmodel2 <- read.csv("Z:/sl2763_Samantha Lau/Projects/2019 FFAR project/Growth curves/R Code/Prelimary_Growth_Curve_Dave_Code_SL_071219.csv")
names(growthmodel2) <- c("t","Isolate","LOG10N","Count")

R100084 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0084", ]

#Run Buchanan growth model
R100084modelb <- nls(buchanan, R100084,
                    list(lag = 6, mumax=1, LOG10N0=3, LOG10Nmax=8))

#Run Gompertz growth model
R100084modelg <- nls(gompertzm, R100084,
                    list(lag = 9, mumax=1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R100084modelba <- nls(baranyi, R100084,
                     list(lag = 9, mumax=1, LOG10N0=3, LOG10Nmax=8))

BIC(R100084modelb)
BIC(R100084modelg)
BIC(R100084modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100084$t, R100084$LOG10N, xlim = c(0,360), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100084$t),max(R100084$t)))
lines(new$t,predict(R100084modelg,newdata=new))

overview(R100084modelg)

R100990 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0990", ]

#Run Buchanan growth model
R100990modelb <- nls(buchanan, R100990,
                     list(lag = 6, mumax=1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100990modelg <- nls(gompertzm, R100990,
                    list(lag = 6, mumax=1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R100990modelba <- nls(baranyi, R100990,
                     list(lag = 6, mumax=1, LOG10N0=3, LOG10Nmax=8))

BIC(R100990modelb)
BIC(R100990modelg)
BIC(R100990modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100990$t, R100990$LOG10N, xlim = c(0,360), ylim = c(3,8), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100990$t),max(R100990$t)))
lines(new$t,predict(R100990modelb,newdata=new))

overview(R100990modelg)

R103286 <- growthmodel2[growthmodel2$Isolate=="FSL R10-3286", ]

#Run Buchanan growth model
R103286modelb <- nls(buchanan, R103286,
                    list(lag = 25, mumax=1, LOG10N0=4, LOG10Nmax=8))

#Run Gompertz growth model
R103286modelg <- nls(gompertzm, R103286,
                    list(lag = 2, mumax=1, LOG10N0=4, LOG10Nmax=8))

#Run Baranyi growth model
R103286modelba <- nls(baranyi, R103286,
                     list(lag = 2, mumax=1, LOG10N0=4, LOG10Nmax=8))

BIC(R103286modelb)
BIC(R103286modelg)
BIC(R103286modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R103286$t, R103286$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R103286$t),max(R103286$t)))
lines(new$t,predict(R103286modelb,newdata=new))

overview(R103286modelg)

A50030 <- growthmodel2[growthmodel2$Isolate=="A5-0030", ]

#Run Buchanan growth model
A50030modelb <- nls(buchanan, A50030,
                    list(lag = 2, mumax=1, LOG10N0=3, LOG10Nmax=9))

#Run Gompertz growth model
A50030modelg <- nls(gompertzm, A50030,
                    list(lag = 2, mumax=1, LOG10N0=3, LOG10Nmax=9))

#Run Baranyi growth model
A50030modelba <- nls(baranyi, A50030,
                     list(lag = 2, mumax=1, LOG10N0=3, LOG10Nmax=9))

BIC(A50030modelb)
BIC(A50030modelg)
BIC(A50030modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Buchanan has lower BIC value
#a rough plot to see how our model looks
plot(A50030$t, A50030$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(A50030$t),max(A50030$t)))
lines(new$t,predict(A50030modelb,newdata=new))

A50030 <- growthmodel2[growthmodel2$Isolate=="A5-0030", ]

#Run Buchanan growth model
A50030modelb <- nls(buchanan, A50030,
                    list(lag = 2, mumax=1, LOG10N0=3, LOG10Nmax=9))

#Run Gompertz growth model
A50030modelg <- nls(gompertzm, A50030,
                    list(lag = 2, mumax=1, LOG10N0=3, LOG10Nmax=9))

#Run Baranyi growth model
A50030modelba <- nls(baranyi, A50030,
                     list(lag = 2, mumax=1, LOG10N0=3, LOG10Nmax=9))

BIC(A50030modelb)
BIC(A50030modelg)
BIC(A50030modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Buchanan has lower BIC value
#a rough plot to see how our model looks
plot(A50030$t, A50030$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(A50030$t),max(A50030$t)))
lines(new$t,predict(A50030modelb,newdata=new))

R70277 <- growthmodel2[growthmodel2$Isolate=="R7-0277", ]

#Run Buchanan growth model
R70277modelb <- nls(buchanan, R70277,
                    list(lag = 19, mumax=1, LOG10N0=5, LOG10Nmax=8))

#Run Gompertz growth model
R70277modelg <- nls(gompertzm, R70277,
                    list(lag = 19, mumax=1, LOG10N0=5, LOG10Nmax=8))

#Run Baranyi growth model
R70277modelba <- nls(baranyi, R70277,
                     list(lag = 19, mumax=1, LOG10N0=5, LOG10Nmax=8))

BIC(R70277modelb)
BIC(R70277modelg)
BIC(R70277modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R70277$t, R70277$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R70277$t),max(R70277$t)))
lines(new$t,predict(R70277modelb,newdata=new))

J30153 <- growthmodel2[growthmodel2$Isolate=="J3-0153", ]

#Run Buchanan growth model
J30153modelb <- nls(buchanan, J30153,
                    list(lag = 1, mumax=1, LOG10N0=3, LOG10Nmax=7))

#Run Gompertz growth model
J30153modelg <- nls(gompertzm, J30153,
                    list(lag = 9, mumax=1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
J30153modelba <- nls(baranyi, J30153,
                     list(lag = 9, mumax=1, LOG10N0=3, LOG10Nmax=8))

BIC(J30153modelb)
BIC(J30153modelg)
BIC(J30153modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(J30153$t, J30153$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(J30153$t),max(J30153$t)))
lines(new$t,predict(J30153modelb,newdata=new))

