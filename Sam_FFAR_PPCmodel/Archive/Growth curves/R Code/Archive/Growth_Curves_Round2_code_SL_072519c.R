#growth curves round 2
#edited by SL on 7/25/19

library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
#insert growth model data
growthmodel2 <- read.csv("Growth_Curve_Round2_allbio_SL_071919.csv") #reading into the
#dataframe
#growthmodel2 <- read.csv("Growth_Curve_Round2_allbio_SL_071919.csv", header=TRUE)
#this makes all the headers that I have in the CSV file as my column names
names(growthmodel2) <- c("t","Isolate","LOG10N","Count", "BioRep") #assigning the names on the file
#even though its the same name

R100531 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0531"& growthmodel2$BioRep=="1", ] #the ending comma
#means to give everything in that row, it finds every run that matches R10-531, IF
#i put log count there, it will give me only the log in that row
#R100531$Round <- "Two"

#Run Buchanan growth model
R100531modelb <- nls(buchanan, R100531,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100531modelg <- nls(gompertzm, R100531,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#non linear system=nls, fitting data, second is data im using, third is starting values

#Run Baranyi growth model
R100531modelba <- nls(baranyi, R100531,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

BICR100531modelb<-BIC(R100531modelb)
BICR100531modelg<-BIC(R100531modelg)
BICR100531modelba<-BIC(R100531modelba)
BICR100531model<- rbind(BICR100531modelb, BICR100531modelg, BICR100531modelba)

#baysian information crtieria and thats calling the information of the model

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100531$t, R100531$LOG10N, xlim = c(0,360), ylim = c(3,9), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100531$t),max(R100531$t)))
lines(new$t,predict(R100531modelg,newdata=new))
lines(new$t,predict(R100531modelb,newdata = new), col="deepskyblue")
#create new data frame and it fills it with the minimum and amximum of the t
#you use this to make the line of the model i have
#input it into the model with the line function, put in the t values (the x values)
#then the y values you predict the line using the model you just built
#spits out values from the model in g

#in plot, first value is the x value, then my value is the log number, then setting 
# limits for the x and y axis, then we want to label the axes
#pch= is the symbol and is a dot
overview(R100531modelb)
R100531bioone<-coef(R100531modelb)

#rbind(R100531bioone, ETC, ETC), THIS WILL COMBINE ALL THE VECTORS TOGETHER
#THEN I CAN WRITE THIS AS A CSV

R101099 <- growthmodel2[growthmodel2$Isolate=="FSL R10-1099"& growthmodel2$BioRep=="1", ]

#Run Buchanan growth model
R101099modelb <- nls(buchanan, R101099,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R101099modelg <- nls(gompertzm, R101099,
                    list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R101099modelba <- nls(baranyi, R101099,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

BICR101099modelb<-BIC(R101099modelb)
BICR101099modelg<-BIC(R1010991modelg)
BICR101099modelba<-BIC(R101099modelba)
BICR101099model<- rbind(BICR101099modelb, BICR101099modelg, BICR101099modelba)

#BUCHANAN IS 3 PHASE LINEAR, SO IT DOES 3 STRAIGHT LINES SO IT WILL PROB BE THE WORSE

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R101099$t, R101099$LOG10N, xlim = c(0,360), ylim = c(3,8), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R101099$t),max(R101099$t)))
lines(new$t,predict(R101099modelb,newdata=new))

overview(R101099modelb)
R101099bioone<-coef(R101099modelb)

R101113 <- growthmodel2[growthmodel2$Isolate=="FSL R10-1113"& growthmodel2$BioRep=="1", ]

#Run Buchanan growth model
R101113modelb <- nls(buchanan, R101113,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R101113modelg <- nls(gompertzm, R101113,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R101113modelba <- nls(baranyi, R101113,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))


BICR101113modelb<-BIC(R101113modelb)
BICR101113modelg<-BIC(R101113modelg)
BICR101113modelba<-BIC(R101113modelba)
BICR101113model<- rbind(BICR101113modelb, BICR101113modelg, BICR101113modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R101113$t, R101113$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R101113$t),max(R101113$t)))
lines(new$t,predict(R101113modelb,newdata=new))

overview(R101113modelg)
R101113bioone<-coef(R101113modelb)

R100941 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0941"& growthmodel2$BioRep=="1", ]

#Run Buchanan growth model
R100941modelb <- nls(buchanan, R100941,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100941modelg <- nls(gompertzm, R100941,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R100941modelba <- nls(baranyi, R100941,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

BIC(R100941modelb)
BIC(R100941modelg)
BIC(R100941modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100941$t, R100941$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100941$t),max(R100941$t)))
lines(new$t,predict(R100941modelb,newdata=new))

overview(R100941modelb)
R100941bioone<-coef(R100941modelb)


R100054 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0054"& growthmodel2$BioRep=="1", ]

#Run Buchanan growth model
R100054modelb <- nls(buchanan, R100054,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100054modelg <- nls(gompertzm, R100054,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R100054modelba <- nls(baranyi, R100054,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

BIC(R100054modelb)
BIC(R100054modelg)
BIC(R100054modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100054$t, R100054$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100054$t),max(R100054$t)))
lines(new$t,predict(R100054modelb,newdata=new))

overview(R100054modelb)
R100054bioone<-coef(R100054modelb)


R100151 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0151"& growthmodel2$BioRep=="1", ]

#Run Buchanan growth model
R100151modelb <- nls(buchanan, R100151,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100151modelg <- nls(gompertzm, R100151,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R100151modelba <- nls(baranyi, R100151,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

BIC(R100151modelb)
BIC(R100151modelg)
BIC(R100151modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100151$t, R100151$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100151$t),max(R100151$t)))
lines(new$t,predict(R100151modelb,newdata=new))

overview(R100151modelb)
R100151bioone<-coef(R100151modelb)



R100553 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0553"& growthmodel2$BioRep=="1", ]

#Run Buchanan growth model
R100553modelb <- nls(buchanan, R100553,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100553modelg <- nls(gompertzm, R100553,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R100553modelba <- nls(baranyi, R100553,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

BICR100553modelb<-BIC(R100553modelb)
BICR100553modelg<-BIC(R100553modelg)
BICR100553modelba<-BIC(R100553modelba)
BICR100553model<- rbind(BICR100553modelb, BICR100553modelg, BICR100553modelba)

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100553$t, R100553$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100553$t),max(R100553$t)))
lines(new$t,predict(R100553modelb,newdata=new))

overview(R100553modelb)
R100553bioone<-coef(R100553modelb)

Growthparametersroundtwo<- rbind(R100531bioone, R101099bioone, R101113bioone, R100941bioone, R100054bioone,  R100553bioone)
write.csv(Growthparametersroundtwo, "Growth_parameters_round_two_072519.csv")
