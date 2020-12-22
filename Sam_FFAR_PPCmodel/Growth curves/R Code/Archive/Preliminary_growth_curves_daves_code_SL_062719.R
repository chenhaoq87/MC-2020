library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
#insert growth model data
growthmodel2 <- read.csv("Prelimary_Growth_Curve_Dave_Code_SL_062719.csv")
names(growthmodel2) <- c("t","Isolate","LOG10N","Count")

R100084 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0084", ]

#Run Buchanan growth model
R100084modelb <- nls(buchanan, R100084,
                    list(lag = 9, mumax=1, LOG10N0=3, LOG10Nmax=8))

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
plot(H80287$t, H80287$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(H80287$t),max(H80287$t)))
lines(new$t,predict(H80287modelb,newdata=new))

J30123 <- growthmodel2[growthmodel2$Isolate=="J3-0123", ]

#Run Buchanan growth model
J30123modelb <- nls(buchanan, J30123,
                    list(lag = 6, mumax=1, LOG10N0=3, LOG10Nmax=8))

#Run Gompertz growth model
J30123modelg <- nls(gompertzm, J30123,
                    list(lag = 6, mumax=1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
J30123modelba <- nls(baranyi, J30123,
                     list(lag = 6, mumax=1, LOG10N0=3, LOG10Nmax=8))

BIC(J30123modelb)
BIC(J30123modelg)
BIC(J30123modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(J30123$t, J30123$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(J30123$t),max(J30123$t)))
lines(new$t,predict(J30123modelb,newdata=new))

H80237 <- growthmodel2[growthmodel2$Isolate=="H8-0237", ]

#Run Buchanan growth model
H80237modelb <- nls(buchanan, H80237,
                    list(lag = 2, mumax=1, LOG10N0=4, LOG10Nmax=8))

#Run Gompertz growth model
H80237modelg <- nls(gompertzm, H80237,
                    list(lag = 2, mumax=1, LOG10N0=4, LOG10Nmax=8))

#Run Baranyi growth model
H80237modelba <- nls(baranyi, H80237,
                     list(lag = 2, mumax=1, LOG10N0=4, LOG10Nmax=8))

BIC(H80237modelb)
BIC(H80237modelg)
BIC(H80237modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(H80237$t, H80237$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(H80237$t),max(H80237$t)))
lines(new$t,predict(H80237modelb,newdata=new))

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

