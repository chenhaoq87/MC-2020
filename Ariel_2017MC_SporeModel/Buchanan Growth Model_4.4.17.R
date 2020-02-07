#Fitting Growth Parameters to Raw Data
#Ariel Buehler, Cornell University, ajb466@cornell.edu
setwd("/Users/arielbuehler/Documents/Documents/Cornell/Research/Sporeformer Characterization/Growth Curve R Work/")
library(nlstools)
library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
library(hydroGOF)#Used for RMSE values
#insert growth model data
growthmodel <- read.csv("growthmodel.csv")
names(growthmodel) <- c("t","Isolate","LOG10N","Count")

J30120 <- growthmodel[growthmodel$Isolate=="J3-0120 ", ]

#Run Buchanan growth model
J30120modelb <- nls(buchanan, J30120,
                   list(lag = 4, mumax=1, LOG10N0=3, LOG10Nmax=8))

#Run Gompertz growth model
J30120modelg <- nls(gompertzm, J30120,
                   list(lag = 4, mumax=1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
J30120modelba <- nls(baranyi, J30120,
                    list(lag = 4, mumax=1, LOG10N0=3, LOG10Nmax=8))

BIC(J30120modelb)
BIC(J30120modelg)
BIC(J30120modelba)

#RMSE calculation
J30120_fittedvalues <- fitted.values(J30120modelb)
rmse(J30120_fittedvalues, J30120$LOG10N)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Gomperts has lower BIC value (so you know: 14.24597 vs. 19.27125)
#a rough plot to see how our model looks
plot(J30120$t, J30120$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(J30120$t),max(J30120$t)))
lines(new$t,predict(J30120modelb,newdata=new))


#just for fun, let's see how the gompertzm looks
lines(new$t,predict(J30120modelg,newdata=new),lty=2)
legend(15, 5, c("Buchanan", "Gompertz"), lty = c(1,2))


#adding baranyi model line to the graph
lines(new$t, predict(J30120modelba, newdata=new), lty=2)
legend(15, 5, c("Buchanan", "Gompertz", "Baranyi"), lty=c(1,2,3))

R50213 <- growthmodel[growthmodel$Isolate=="R5-0213 ", ]

#Run Buchanan growth model
R50213modelb <- nls(buchanan, R50213,
                    list(lag = 4, mumax=1, LOG10N0=3, LOG10Nmax=8))

R50213modelbar <- nls(baranyi, R50213, 
                      list(lag = 4, mumax=1, LOG10N0 =3, LOG10Nmax=8))

R50213modelg <- nls(gompertzm, R50213,
                    list(lag=4, mumax=1, LOG10N0=3, LOG10Nmax=8))

BIC(R50213modelb)
BIC(R50213modelbar)
BIC(R50213modelg)

plot(R50213$t, R50213$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R50213$t),max(R50213$t)))
lines(new$t,predict(R50213modelb,newdata=new))

H70687 <- growthmodel[growthmodel$Isolate=="H7-0687 ", ]

#Run Buchanan growth model
H70687modelb <- nls(buchanan, H70687,
                    list(lag = 4, mumax=1, LOG10N0=3, LOG10Nmax=8))

H70687modelbar <- nls(baranyi, H70687,
                    list(lag = 4, mumax=1, LOG10N0=3, LOG10Nmax=8))

H70687modelg <- nls(gompertzm, H70687,
                    list(lag = 4, mumax=1, LOG10N0=3, LOG10Nmax=8))

BIC(H70687modelb)
BIC(H70687modelbar)
BIC(H70687modelg)

plot(H70687$t, H70687$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(H70687$t),max(H70687$t)))
lines(new$t,predict(H70687modelb,newdata=new))

W80169 <- growthmodel[growthmodel$Isolate=="W8-0169 ", ]

#Run Buchanan growth model
W80169modelb <- nls(buchanan, W80169,
                    list(lag = 12, mumax=1, LOG10N0=4, LOG10Nmax=8))

W80169modelbar <- nls(baranyi, W80169,
                    list(lag = 12, mumax=1, LOG10N0=4, LOG10Nmax=8))

W80169modelg <- nls(gompertzm, W80169,
                    list(lag = 12, mumax=1, LOG10N0=4, LOG10Nmax=8))

BIC(W80169modelb)
BIC(W80169modelbar)
BIC(W80169modelg)

plot(W80169$t, W80169$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(W80169$t),max(W80169$t)))
lines(new$t,predict(W80169modelb,newdata=new))

##Trial two growth models
growthmodel2 <- read.csv("growthmodeltrial2_nooutlier_4.4.17.csv")
names(growthmodel2) <- c("t","Isolate","LOG10N","Count")

H80287 <- growthmodel2[growthmodel2$Isolate=="H8-0287", ]

#Run Buchanan growth model
H80287modelb <- nls(buchanan, H80287,
                    list(lag = 9, mumax=1, LOG10N0=3, LOG10Nmax=8))

#Run Gompertz growth model
H80287modelg <- nls(gompertzm, H80287,
                    list(lag = 9, mumax=1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
H80287modelba <- nls(baranyi, H80287,
                     list(lag = 9, mumax=1, LOG10N0=3, LOG10Nmax=8))

BIC(H80287modelb)
BIC(H80287modelg)
BIC(H80287modelba)

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
                    list(lag = 3, mumax=1, LOG10N0=3.5, LOG10Nmax=6.5))

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
plot(A50030$t, A50030$LOG10N, xlim = c(0,27), ylim = c(3,8), xlab = "Time (days)", ylab = expression(paste(LOG[10],CFU/mL)), pch = 20)
new = data.frame(t = seq(min(A50030$t),max(A50030$t)))
lines(new$t,predict(A50030modelg,newdata=new), lty=c(1))
legend(15, 5, c("Gompertz"), lty=c(1))

#just for fun, let's see how the gompertzm looks
lines(new$t,predict(A50030modelba,newdata=new),lty=3, col=34, lwd=c(2))
legend(15, 5, c("Gompertz", "Baranyi"), lty = c(1,3), lwd=c(2))


#adding baranyi model line to the graph
lines(new$t, predict(A50030modelb, newdata=new), lty=4, col=28, lwd=c(2))
legend(15, 5, c("Gompertz", "Baranyi", "Buchanan"), lty=c(1,3,2), col=c(153, 34, 28), lwd=c(2))

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
AIC(A50030modelb)
AIC(A50030modelg)
AIC(A50030modelba)

#RMSE calculation
A50030modelb_fittedvalues <- fitted.values(A50030modelb)
rmse(A50030modelb_fittedvalues, A50030$LOG10N)
A50030modelg_fittedvalues <- fitted.values(A50030modelg)
rmse(A50030modelg_fittedvalues, A50030$LOG10N)
A50030modelba_fittedvalues <- fitted.values(A50030modelba)
rmse(A50030modelba_fittedvalues, A50030$LOG10N)

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








#summary plot of all Buchanan models fittedYAY!!!!!! FINALLLLYYYYYYYYY!!!!!
plot(new$t, predict(A50030modelb, newdata=new), type= "l", xlab = "Time (days)", ylab=expression(paste(LOG[10],CFU/mL)), pch=20, lwd=1.5)
lines(new$t, predict(W80169modelb, newdata=new), col="red", lwd=1.5)
lines(new$t, predict(J30120modelb, newdata=new), col="blue", lwd=1.5)
lines(new$t, predict(H70687modelb, newdata=new), col="green4", lwd=1.5)
lines(new$t, predict(H80287modelb, newdata=new), col="blueviolet", lwd=1.5)
lines(new$t, predict(J30123modelb, newdata=new), col="cyan", lwd=1.5)
lines(new$t, predict(H80237modelb, newdata=new), col="darkorange", lwd=1.5)
lines(new$t, predict(R50213modelb, newdata=new), col="sienna", lwd=1.5)
lines(new$t, predict(R70277modelb, newdata=new), col="deeppink", lwd=1.5)
legend(22, 5.5, c("AT179", "AT61", "AT340", "AT3", "AT100", "AT513", "AT15", "AT17", "AT45"), 
       lty=c(1,1,1,1), lwd=c(1.5, 1.5, 1.5, 1.5), 
       col=c("black", "red", "blue", "green4", "blueviolet", "cyan", "darkorange", "sienna", "deeppink"), cex=.5)
