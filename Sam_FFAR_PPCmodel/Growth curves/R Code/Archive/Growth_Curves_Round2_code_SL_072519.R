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

BIC(R100531modelb)
BIC(R100531modelg)
BIC(R100531modelba)

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

BIC(R101099modelb)
BIC(R101099modelg)
BIC(R101099modelba)

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

BIC(R101113modelb)
BIC(R101113modelg)
BIC(R101113modelba)

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
