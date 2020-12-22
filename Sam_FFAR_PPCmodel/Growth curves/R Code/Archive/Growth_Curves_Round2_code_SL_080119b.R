#growth curves round 2
#edited by SL on 7/30/19

library(nlsMicrobio) #Used to run microbial growth models
library(nlme) #Used to get BIC values
#insert growth model data
growthmodel2 <- read.csv("Growth_Curve_Round2_allbio_SL_071919.csv", stringsAsFactors = F) #reading into the
#dataframe
#growthmodel2 <- read.csv("Growth_Curve_Round2_allbio_SL_071919.csv", header=TRUE)
#this makes all the headers that I have in the CSV file as my column names
names(growthmodel2) <- c("t","Isolate","LOG10N","Count", "BioRep") #assigning the names on the file
#even though its the same name

for (i in seq(1,3)) { 

R100531 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0531"& growthmodel2$BioRep==i, ] #the ending comma
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
overview(R100531modelg)
R100531b<-coef(R100531modelb)

R100531b["Isolate"] <- R100531$Isolate[1]
R100531b["BIC_b"]<-BIC(R100531modelb)
R100531b["BIC_g"]<-BIC(R100531modelg)
R100531b["BIC_ba"]<-BIC(R100531modelba)

R100531g<-coef(R100531modelg)

R100531g["Isolate"] <- R100531$Isolate[1]
R100531g["BIC_b"]<-BIC(R100531modelb)
R100531g["BIC_g"]<-BIC(R100531modelg)
R100531g["BIC_ba"]<-BIC(R100531modelba)

R100531ba<-coef(R100531modelba)

R100531ba["Isolate"] <- R100531$Isolate[1]
R100531ba["BIC_b"]<-BIC(R100531modelb)
R100531ba["BIC_g"]<-BIC(R100531modelg)
R100531ba["BIC_ba"]<-BIC(R100531modelba)

R101099 <- growthmodel2[growthmodel2$Isolate=="FSL R10-1099"& growthmodel2$BioRep==i, ]

#Run Buchanan growth model
R101099modelb <- nls(buchanan, R101099,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R101099modelg <- nls(gompertzm, R101099,
                    list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))

#Run Baranyi growth model
R101099modelba <- nls(baranyi, R101099,
                     list(lag = 6, mumax=0.1, LOG10N0=3, LOG10Nmax=8))



#BUCHANAN IS 3 PHASE LINEAR, SO IT DOES 3 STRAIGHT LINES SO IT WILL PROB BE THE WORSE

#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R101099$t, R101099$LOG10N, xlim = c(0,360), ylim = c(3,8), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R101099$t),max(R101099$t)))
lines(new$t,predict(R101099modelb,newdata=new))

overview(R101099modelb)
R101099b<-coef(R101099modelb)

R101099b["Isolate"] <- R101099$Isolate[1]
R101099b["BIC_b"]<-BIC(R100531modelb)
R101099b["BIC_g"]<-BIC(R100531modelg)
R101099b["BIC_ba"]<-BIC(R100531modelba)

R101099g<-coef(R101099modelg)

R101099g["Isolate"] <- R101099$Isolate[1]
R101099g["BIC_b"]<-BIC(R100531modelb)
R101099g["BIC_g"]<-BIC(R100531modelg)
R101099g["BIC_ba"]<-BIC(R100531modelba)

R101099ba<-coef(R101099modelba)

R101099ba["Isolate"] <- R101099$Isolate[1]
R101099ba["BIC_b"]<-BIC(R100531modelb)
R101099ba["BIC_g"]<-BIC(R100531modelg)
R101099ba["BIC_ba"]<-BIC(R100531modelba)

R101113 <- growthmodel2[growthmodel2$Isolate=="FSL R10-1113"& growthmodel2$BioRep==i, ]

#Run Buchanan growth model
R101113modelb <- nls(buchanan, R101113,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R101113modelg <- nls(gompertzm, R101113,
                    list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R101113modelba <- nls(baranyi, R101113,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))



#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R101113$t, R101113$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R101113$t),max(R101113$t)))
lines(new$t,predict(R101113modelb,newdata=new))

overview(R101113modelg)
R101113b<-coef(R101113modelb)

R101113b["Isolate"] <- R101113$Isolate[1]
R101113b["BIC_b"]<-BIC(R101113modelb)
R101113b["BIC_g"]<-BIC(R101113modelg)
R101113b["BIC_ba"]<-BIC(R101113modelba)

R101113g<-coef(R101113modelb)

R101113g["Isolate"] <- R101113$Isolate[1]
R101113g["BIC_b"]<-BIC(R101113modelb)
R101113g["BIC_g"]<-BIC(R101113modelg)
R101113g["BIC_ba"]<-BIC(R101113modelba)

R101113ba<-coef(R101113modelba)

R101113ba["Isolate"] <- R101113$Isolate[1]
R101113ba["BIC_b"]<-BIC(R101113modelb)
R101113ba["BIC_g"]<-BIC(R101113modelg)
R101113ba["BIC_ba"]<-BIC(R101113modelba)



R100941 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0941"& growthmodel2$BioRep==i, ]

#Run Buchanan growth model
R100941modelb <- nls(buchanan, R100941,
                     list(lag = 1, mumax=0.07, LOG10N0=3, LOG10Nmax=8), 
                     control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))

#Run Gompertz growth model
R100941modelg <- nls(gompertzm, R100941,
                     list(lag = 8, mumax=0.07, LOG10N0=3, LOG10Nmax=8), 
                     control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))

#Run Baranyi growth model
R100941modelba <- nls(baranyi, R100941,
                      list(lag = 8, mumax=0.07, LOG10N0=3, LOG10Nmax=8), 
                      control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))



#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100941$t, R100941$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100941$t),max(R100941$t)))
lines(new$t,predict(R100941modelb,newdata=new))

overview(R100941modelb)
R100941b<-coef(R100941modelb)

R100941b["Isolate"] <- R100941$Isolate[1]
R100941b["BIC_b"]<-BIC(R100941modelb)
R100941b["BIC_g"]<-BIC(R100941modelg)
R100941b["BIC_ba"]<-BIC(R100941modelba)

R100941g<-coef(R100941modelg)

R100941g["Isolate"] <- R100941$Isolate[1]
R100941g["BIC_b"]<-BIC(R100941modelb)
R100941g["BIC_g"]<-BIC(R100941modelg)
R100941g["BIC_ba"]<-BIC(R100941modelba)

R100941ba<-coef(R100941modelba)

R100941ba["Isolate"] <- R100941$Isolate[1]
R100941ba["BIC_b"]<-BIC(R100941modelb)
R100941ba["BIC_g"]<-BIC(R100941modelg)
R100941ba["BIC_ba"]<-BIC(R100941modelba)

R100054 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0054"& growthmodel2$BioRep==i, ]

#Run Buchanan growth model
R100054modelb <- nls(buchanan, R100054,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100054modelg <- nls(gompertzm, R100054,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R100054modelba <- nls(baranyi, R100054,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))


#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100054$t, R100054$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100054$t),max(R100054$t)))
lines(new$t,predict(R100054modelb,newdata=new))

overview(R100054modelg)
R100054b<-coef(R100054modelb)

R100054b["Isolate"] <- R100054$Isolate[1]
R100054b["BIC_b"]<-BIC(R100054modelb)
R100054b["BIC_g"]<-BIC(R100054modelg)
R100054b["BIC_ba"]<-BIC(R100054modelba)

R100054g<-coef(R100054modelg)

R100054g["Isolate"] <- R100054$Isolate[1]
R100054g["BIC_b"]<-BIC(R100054modelb)
R100054g["BIC_g"]<-BIC(R100054modelg)
R100054g["BIC_ba"]<-BIC(R100054modelba)

R100054ba<-coef(R100054modelba)

R100054ba["Isolate"] <- R100054$Isolate[1]
R100054ba["BIC_b"]<-BIC(R100054modelb)
R100054ba["BIC_g"]<-BIC(R100054modelg)
R100054ba["BIC_ba"]<-BIC(R100054modelba)

R100151 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0151"& growthmodel2$BioRep==i, ]

#Run Buchanan growth model
R100151modelb <- nls(buchanan, R100151,
                     list(lag = 8, mumax=0.07, LOG10N0=3, LOG10Nmax=8), 
                     control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))

#Run Gompertz growth model
R100151modelg <- nls(gompertzm, R100151,
                     list(lag = 8, mumax=0.07, LOG10N0=3, LOG10Nmax=8), 
                     control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))

#Run Baranyi growth model
R100151modelba <- nls(baranyi, R100151,
                      list(lag = 8, mumax=0.07, LOG10N0=3, LOG10Nmax=8), 
                      control = nls.control(maxiter = 100, minFactor = 1/4096, warnOnly = T))



#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100151$t, R100151$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100151$t),max(R100151$t)))
lines(new$t,predict(R100151modelb,newdata=new))

overview(R100151modelb)
R100151b<-coef(R100151modelb)

R100151b["Isolate"] <- R100151$Isolate[1]
R100151b["BIC_b"]<-BIC(R100151modelb)
R100151b["BIC_g"]<-BIC(R100151modelg)
R100151b["BIC_ba"]<-BIC(R100151modelba)

R100151g<-coef(R100151modelg)

R100151g["Isolate"] <- R100151$Isolate[1]
R100151g["BIC_b"]<-BIC(R100151modelb)
R100151g["BIC_g"]<-BIC(R100151modelg)
R100151g["BIC_ba"]<-BIC(R100151modelba)

R100151ba<-coef(R100151modelba)

R100151ba["Isolate"] <- R100151$Isolate[1]
R100151ba["BIC_b"]<-BIC(R100151modelb)
R100151ba["BIC_g"]<-BIC(R100151modelg)
R100151ba["BIC_ba"]<-BIC(R100151modelba)

R100553 <- growthmodel2[growthmodel2$Isolate=="FSL R10-0553"& growthmodel2$BioRep==i, ]

#Run Buchanan growth model
R100553modelb <- nls(buchanan, R100553,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Gompertz growth model
R100553modelg <- nls(gompertzm, R100553,
                     list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))

#Run Baranyi growth model
R100553modelba <- nls(baranyi, R100553,
                      list(lag = 6, mumax=0.1, LOG10N0=1, LOG10Nmax=8))


#These models look about the same but Buchanan is easier to interpret. Use Buchanan.
#Baranyi has lower BIC value
#a rough plot to see how our model looks
plot(R100553$t, R100553$LOG10N, xlim = c(0,360), ylim = c(3,10), xlab = "Time (hours)", ylab = "LOG10N", pch = 20)
new = data.frame(t = seq(min(R100553$t),max(R100553$t)))
lines(new$t,predict(R100553modelb,newdata=new))

overview(R100553modelb)
R100553b<-coef(R100553modelb)
R100553b["Isolate"] <- R100553$Isolate[1]
R100553b["BIC_b"]<-BIC(R100553modelb)
R100553b["BIC_g"]<-BIC(R100553modelg)
R100553b["BIC_ba"]<-BIC(R100553modelba)

R100553g<-coef(R100553modelg)
R100553g["Isolate"] <- R100553$Isolate[1]
R100553g["BIC_b"]<-BIC(R100553modelb)
R100553g["BIC_g"]<-BIC(R100553modelg)
R100553g["BIC_ba"]<-BIC(R100553modelba)

R100553ba<-coef(R100553modelba)
R100553ba["Isolate"] <- R100553$Isolate[1]
R100553ba["BIC_b"]<-BIC(R100553modelb)
R100553ba["BIC_g"]<-BIC(R100553modelg)
R100553ba["BIC_ba"]<-BIC(R100553modelba)

Growthparametersround2b<-as.data.frame(rbind(R100531b, R101099b, R101113b, R100941b, R100054b,  R100151b, R100553b))
Growthparametersround2b$Bio_Replicate <- i
Growthparametersround2b<-Growthparametersround2b[, c(5, 9, 1, 2, 3, 4)]
#colnames(Growthparameteround2b)[3]<-"lagb"
#colnames(Growthparameteround2b)[3]<-paste("lagb")
colnames(Growthparametersround2b)[colnames(Growthparametersround2b)=="lag"] <- "lag_b"
colnames(Growthparametersround2b)[colnames(Growthparametersround2b)=="mumax"] <- "mumax_b"
colnames(Growthparametersround2b)[colnames(Growthparametersround2b)=="LOG10N0"] <- "log10N0_b"
colnames(Growthparametersround2b)[colnames(Growthparametersround2b)=="LOG10Nmax"] <- "log10Nmax_b"

Growthparametersround2ba<-as.data.frame(rbind(R100531ba, R101099ba, R101113ba, R100941ba, R100054ba,  R100151ba, R100553ba))
Growthparametersround2ba$Bio_Replicate <- i
Growthparametersround2ba<-Growthparametersround2ba[, c(1, 2, 3, 4)]
colnames(Growthparametersround2ba)[colnames(Growthparametersround2ba)=="lag"] <- "lag_ba"
colnames(Growthparametersround2ba)[colnames(Growthparametersround2ba)=="mumax"] <- "mumax_ba"
colnames(Growthparametersround2ba)[colnames(Growthparametersround2ba)=="LOG10N0"] <- "log10N0_ba"
colnames(Growthparametersround2ba)[colnames(Growthparametersround2ba)=="LOG10Nmax"] <- "log10Nmax_ba"

Growthparametersround2g<-as.data.frame(rbind(R100531g, R101099g, R101113g, R100941g, R100054g,  R100151g, R100553g))
Growthparametersround2g$Bio_Replicate <- i
Growthparametersround2g<-Growthparametersround2g[, c(1, 2, 3, 4, 6, 7, 8)]
colnames(Growthparametersround2g)[colnames(Growthparametersround2g)=="lag"] <- "lag_g"
colnames(Growthparametersround2g)[colnames(Growthparametersround2g)=="mumax"] <- "mumax_g"
colnames(Growthparametersround2g)[colnames(Growthparametersround2g)=="LOG10N0"] <- "log10N0_g"
colnames(Growthparametersround2g)[colnames(Growthparametersround2g)=="LOG10Nmax"] <- "log10Nmax_g"


#Growthparametersall<-as.data.frame(rbind(Growthparametersround2b, Growthparametersround2ba, Growthparametersround2g))

Growthparametersall<-as.data.frame(cbind(Growthparametersround2b, Growthparametersround2ba, Growthparametersround2g))
#colnames(Growthparametersall)[7]<="lag_ba"

if (i==1) {
  write.table(Growthparametersall, "Growthparametersall_080119.csv", row.names = FALSE, sep=",",append=FALSE)
  #instead of rewriting the file each time, it adds things to the end
  #sep tells it what to write in between , a comma means different column
  
  
}
else {

write.table(Growthparametersall, "Growthparametersall_080119.csv", sep=",",  row.names = FALSE, col.names=FALSE, append=TRUE)
  
#instead of rewriting the file each time, it adds things to the end
#sep tells it what to write in between , a comma means different column
}
}
#Growthparametersroundtwo<- rbind(R100531bioone, R101099bioone, R101113bioone, R100941bioone, R100054bioone,  R100553bioone)
#write.csv(Growthparametersroundtwo, "Growth_parameters_round_two_072519.csv")
