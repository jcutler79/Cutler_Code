
# Note: I have not cleaned up most of this code--a lot of it is still
## kind of rambling or jumbled, with fits and starts, with no clear, annotated
## flow or logic. Sorry about that! If you have any questions let me know.

# Survival Analysis


# Online and downloaded sources
## Textbook
## drizopoulos.com "EP03: Survival Analysis in R Companion" - http://www.drizopoulos.com/courses/emc/ep03_%20survival%20analysis%20in%20r%20companion#estimators-of-the-survival-function
## For some reason the Nelson-Aalen gets mentioned as something that can't be calculated DIRECTLY in the survival package: https://folk.ntnu.no/bo/TMA4275/Download/R.tutorialDiez.pdf
## For non-cumulative hazard plots: https://rdrr.io/cran/rms/man/survplot.html (rms package)
## Finalfit: Survival analysis with strata, clusters, frailties and competing risks in in Finalfit - https://www.r-bloggers.com/survival-analysis-with-strata-clusters-frailties-and-competing-risks-in-in-finalfit/
## survplot, from rms - https://rdrr.io/cran/rms/man/survplot.html
## log hazard ratio
## Drawing survival curves using ggplot2: https://rpkgs.datanovia.com/survminer/reference/ggsurvplot.html






# Libraries
library(ggplot2)
library(data.table)
library(tidyr)
library(broom)

# library(grDevices)
library(gridExtra)
library(wrapr)


# Must have survival packages
library(survival)
library(KMsurv)         # What is this for? one source on it: http://biostat.mc.vanderbilt.edu/wiki/pub/Main/QingxiaChen/Lab12.txt
library(survminer)      # For ggsurvplot
library(survMisc)       # for log rank trend test (use comp on a ten object) - screws with autoplot though - hey it even warns that that's what it does!
library(ranger)         # install.packages("ranger")     - What are these two for?
library(ggfortify)      # install.packages("ggfortify")  - What are these two for?
library(rms)            # Regression Modeling Strategies - includes survplot for (non-cumulative, i.e. instantaneous) hazard plots!
## https://rdrr.io/cran/rms/man/survplot.html


library(muhaz)          # install.packages("muhaz")
library(discSurv)       # install.packages("discSurv") - For making a life table!
# library(mixPHM)       # This one is too weird to use


# From drizopoulos
library(splines)
library(lattice)
library(JM)


# library(rgl)
# library(kernhaz)
# library(survPresmooth)







# Grades
# HW           35%
# Midterm      25%
# Presentation 15%
# Project      25%
hw = .88
mid = 56/60
pres = .9
proj = .9
hw.w = 35
mid.w = 25
pres.w = 15
proj.w = 25
sum(hw*hw.w,
    mid*mid.w,
    pres*pres.w,
    proj*proj.w)






# Datasets in my Survival Desktop folder
## 6-mp                                               = leukemia.csv
## Time to discontinuation of IUD                     =
## Survival of multiple myeloma                       = multiple_myeloma.csv
## Prognosis for women with breast cancer             =
## Comparison of two treatments for prostate cancer   = Prostate_Cancer.csv
## Recurrence of bladder cancer                       = bladder_cancer.csv
## Survival times of patients with melanoma           = melanoma.csv
## PBC  (primary billiary cirrhosis)                  = pbcirrhosis.csv          










# My custom functions

# Now I can do my own freaking LR test to compare nested models!
## This was only necessary once I tried to compare a flexsurvreg gamma model
## to any of its nested models, lognormal, weibull, and exponential. anova()
## doesn't know how to work with a flexsurvreg object. Thank goodness it's 
## a simple fix!
myLRT = function(loglik1,loglik2, df1,df2){ # df1 and df2 refer to the DEGREES OF FREEDOM; not to dataframes
  neg2logL1 = -2*loglik1
  neg2logL2 = -2*loglik2
  lrt = abs(neg2logL1-neg2logL2)
  finaldf = abs(df1-df2)
  pval = 1-pchisq(lrt,finaldf)
  print("LRT ... DF ... P-value")
  return(c(lrt,finaldf,pval))
}


# From survival_week4hw_Cutler
# Function for creating a "log-rank" table without wasting huge 
## amounts of time and increasing risk of errors from doing it by hand:
myLogRankTable = function(groups.1_2.TimeStatus.DF){
  dt = data.table::data.table(groups.1_2.TimeStatus.DF)
  # 1) order by time
  dtfull = dt[order(time)]
  # 2) just events, no censored
  dt1 = dtfull[status == 1]
  # 3) eliminate time duplicates
  dt2 = dt1[,sum(status),by="time"]
  colnames(dt2) = c("time","deaths")
  # don't forget to match the group column when adding it:
  dt2$group = dt1$group[match(dt2$time,dt1$time)]
  dt2 = dt2[,c(2,3,1)]
  # 4) add d1j
  dt2$d1j = NA
  deaths1 = dt2[group == 1,which=TRUE]
  dt2$d1j[deaths1] = dt2$deaths[deaths1]
  dt2$d1j[-deaths1] = 0
  # 5) add n1j
  ## add group1-only censored column to full df:
  cens1 = dtfull[group == 1 & status == 0,which=TRUE]
  dtfull$n1j.cens = NA
  dtfull$n1j.cens[cens1] = 1
  dtfull$n1j.cens[-cens1] = 0
  ## d1j for full df:
  dtfull$d1j = NA
  deaths1 = dtfull[group == 1 & status == 1,which=TRUE]
  dtfull$d1j[deaths1] = 1
  dtfull$d1j[-deaths1] = 0
  ## n1j for full df:
  n1j = vector()
  for (i in 1:nrow(dtfull)){
    n1j[i] = dtfull[group == 1,.N] -
      ( sum(dtfull$d1j[1:i-1]) + sum(dtfull$n1j.cens[1:(i-1)]) )
  }
  dtfull$n1j = n1j
  ## and now match-add an n1j column to dt2:
  dt2$n1j =
    dtfull$n1j[match(dt2$time,dtfull$time)]
  # 6) add d2j
  dt2$d2j = NA
  deaths2 = dt2[group == 2,which=TRUE]
  dt2$d2j[deaths2] = dt2$deaths[deaths2]
  dt2$d2j[-deaths2] = 0
  # 7) add n2j
  ## add group2-only censored column to full df:
  cens2 = dtfull[group == 2 & status == 0,which=TRUE]
  dtfull$n2j.cens = NA
  dtfull$n2j.cens[cens2] = 1
  dtfull$n2j.cens[-cens2] = 0
  ## d2j for full df:
  dtfull$d2j = NA
  deaths2 = dtfull[group == 2 & status == 1,which=TRUE]
  dtfull$d2j[deaths2] = 1
  dtfull$d2j[-deaths2] = 0
  ## n2j for full df:
  n2j = vector()
  for (i in 1:nrow(dtfull)){
    n2j[i] = dtfull[group == 2,.N] -
      ( sum(dtfull$d2j[1:i-1]) + sum(dtfull$n2j.cens[1:(i-1)]) )
  }
  dtfull$n2j = n2j
  ## and now match-add an n2j column to dt2:
  dt2$n2j =
    dtfull$n2j[match(dt2$time,dtfull$time)]
  # 8) add d1j and d2j to get dj
  dt2$dj = dt2$d1j + dt2$d2j
  # 9) add n1j and n2j to get nj
  dt2$nj = dt2$n1j + dt2$n2j
  # 10) get e1j
  dt2$e1j = dt2$n1j*dt2$dj/dt2$nj
  # 11) get v1j
  dt2$v1j = ( dt2$n1j*dt2$n2j*dt2$dj*(dt2$nj-dt2$dj) ) /
    ( (dt2$nj^2)*(dt2$nj - 1) )
  # 12) get O1j - Ej column
  dt2$OminusE = dt2$d1j - dt2$e1j
  X2L = (sum(dt2$OminusE)^2)/sum(dt2$v1j)
  X2L.pvalue =
    1-pchisq(X2L,df=1)
  XL = sum(dt2$OminusE)/sqrt(sum(dt2$v1j))
  print(
    sprintf(
      "The chi-square (X^2 df=L) statistic and p-value are %s and %s, respectively, and XL is %s",
      round(X2L,3),round(X2L.pvalue,3),round(XL,3) ))
  return(list(dt2,dtfull,X2L,X2L.pvalue,XL))
}


## Autoplot a life table survival curve 
## (Using the original dataframe, the life table object, and the interval borders)
myLTautoplot = function(disc.life,intBorders){
  new.disc.df = data.frame(time = intBorders,
                           survival = disc.life$Output$S)
  ggplot(new.disc.df, aes(time,survival)) +
    geom_line() +
    geom_point() +
    coord_cartesian(xlim = c(0,max(intBorders)),ylim = c(0,1)) +
    labs(title = "Life Table Survival Curve",
         x="Time",y="Survival Probability")
  #return(new.disc.df)
}
# plot(l.limits,l.life$Output$S)

## A function for getting to the point where you can use data.table querying on 
## your summary object variables of interest.
summ.to.df = function(summ.object,cols = c(2,3,4,6,9,10)){ 
  mycols = names(summ.object)[cols]
  mydf = 
    lapply(mycols,function(i){
      thedf = data.frame(summ.object[[i]])
      return(thedf)
    })
  mydf = do.call("cbind",mydf)
  colnames(mydf) = mycols
  return(mydf)
}


## Convert quartile survival times and CIs to a dataframe 
quant.to.df = function(quant.object){
  df = rbind(
    quant.object[[1]],
    quant.object[[2]],
    quant.object[[3]]
  )
  df = t(df)
  colnames(df) = c("point estimate","lower","upper")
  df = as.data.frame(df)
  return(df)
}


## Convert a string of a DAT file data to a data frame!
## This is actually really important because you can't even use excel to open 
## a DAT! It messes it up!
myDATstring.to.df = function(mystring,mystringNames){
  mystring = unlist(strsplit(mystring,"\n"))
  mystring = gsub("^\\s+","",mystring)
  mystring = gsub("\\s+"," ",mystring)
  mystring = unlist(strsplit(mystring," "))
  mystringNames = gsub("\\s+"," ",mystringNames)
  mystringNames = unlist(strsplit(mystringNames," "))
  myDF =
    lapply(mystringNames,function(i){
      i = mystring[seq(which(mystringNames == i),
                       length(mystring),
                       length(mystringNames))] 
    })
  myDF = do.call("cbind",myDF)
  myDF = as.data.frame(myDF)
  colnames(myDF) = mystringNames
  for (i in 1:ncol(myDF)){
    myDF[,i] = as.numeric(as.character(myDF[,i]))
  } 
  return(myDF)
}



## Function for creating my own intervalBorders character vector
myintBorders = function(intervalLimits){
  myints = vector()
  j = 0
  for (i in 1:length(intervalLimits)){
    myints[i] = paste0("[",j,",",intervalLimits[i],")")
    j = intervalLimits[i]
  }
  return(myints)
}
myintBorders(c(6,12,24,36))












######################################################################
######################################################################

# First in class exercise:

# write.csv(leuk, file = "/Users/jamescutler/Desktop/Survival/leukemia.csv")
load("/Users/jamescutler/Desktop/Survival/leukemia.csv")

# Make a survfit object:
km.fit = survfit(Surv(time, delta) ~ group, data = leuk)
summary(km.fit)

# Plot it:
autoplot(km.fit) +
  theme_bw()

# Plot the cumulative hazard function:
ggsurvplot(km.fit, conf.int = TRUE,
           risk.table = TRUE, risk.table.col = "strata",
           fun = "cumhaz")


# Fake data:
class = 
  build_frame(
    "time","delta"|
      1   ,1 |
      2   ,0 |
      3   ,1 |
      4   ,1 |
      5   ,1 |
      7   ,0
  )

c.fit = survfit(Surv(time,delta) ~ 1, data = class)
summary(c.fit)
autoplot(c.fit)







######################################################################
######################################################################

# Homework 1

# (See week2hw_survival.R)





######################################################################
######################################################################

# Multiple myeloma DAT file data

# Used my custom myDATstring.to.df function to convert the myeloma DAT file 
## string to a data frame.

# write.csv(myeloma,file = "/Users/jamescutler/Desktop/Survival/multiple_myeloma.csv")
myeloma = read.csv("/Users/jamescutler/Desktop/Survival/multiple_myeloma.csv")

# Make it easy to index by protein just for fun:
mmdt = data.table(myeloma)
mymm = mmdt[protein == 1]

# How you get the numbers in the order you need to make the life table by hand:
mymm[order(time)]     

# Start the analysis with creating the survival object:
mm.fit = survfit(Surv(time,status) ~ 1, data = mymm)
mm.summ = summary(mm.fit)
summ.to.df(mm.summ)

# An autoplot:
autoplot(mm.fit) +
  theme_bw()

# A cumulative hazard function plot:
ggsurvplot(mm.fit, conf.int = TRUE,
           risk.table = TRUE, risk.table.col = "strata",
           fun = "cumhaz")





####################

# Returning now to leukemia:

leuk = read.csv("/Users/jamescutler/Desktop/Survival/leukemia.csv")
leuk = leuk[,-1]
head(leuk)

leuk.fit = survfit(Surv(time,delta) ~ group, data = leuk)
summary(leuk.fit)

autoplot(leuk.fit) +
  theme_bw() +
  labs(title = "Leukemia survival: treatment vs control groups",
       x="Time in months",y="Survival probability")

# This looks just like the Nelson-Aalen estimate plot from the SAS output
## in wk3 lecture from class (slide 20) -- YAY!!!!!!!
ggsurvplot(leuk.fit, 
           conf.int = TRUE, 
           risk.table = TRUE,
           risk.table.col = "strata",     # Color-codes the numbers in the risk table. That's it.
           fun = "cumhaz")

# Now I just need to do two things:
## 1) Make a life table with start age and end age in two adjacent columns
## 2) Make an Epanechnikov kernel-smoothed cumulative hazard function plot

# SAS life table:
build_frame(
  "[lower,","upper)","# failed","# censored","effective n","P(F|?)","C.P. SE","S","failure","S. SE","median residual lifetime"
)

# discSurv life table:
build_frame(
  "[lower, upper)","n","events","dropouts","atRisk","hazard","H. SE","S","S. SE","cumHazard","C.H. SE","margProb = P(F)?"
)

lifeTable(dataSet = leuk, 
          timeColumn = "time",
          censColumn = "delta")


# Try the life table out with the fake dataset:
dummy = 
  build_frame(
    "j" , "c" |
      1 , 1 |   #
      2 , 0 |   ##
      3 , 1 |   ###
      4 , 1 |   ####
      5 , 1 |   #####
    7 , 0     #####
  )

# Discretize the intervals:
intBorders = c(1,3,4,5,8)
dummy.disc = 
  contToDisc(dataSet = dummy, 
             timeColumn = "j", 
             intervalLimits = intBorders)

# Make the life table:
d.disc.life = 
  lifeTable(dummy.disc,
            timeColumn = "timeDisc",
            censColumn = "c")
d.disc.life

# Ha! This is so simple! Just plot a line graph of the life table survival column duh!!

# Make a function to autoplot the life table survival curve:
## (See custom functions at top of notepad)
myLTautoplot(dummy.disc,d.disc.life,intBorders = intBorders)







######################################################################
######################################################################

# Week 3 Homework

# m.names = "age treatment time status"
# melan = 
#   "1 1 19 1
# 1 1 24 0
# 1 1 8 1
# 1 1 17 0
# 1 1 17 0
# 1 1 34 0
# 1 2 27 0
# 1 2 21 0
# 1 2 18 0
# 1 2 16 0
# 1 2 7 1
# 1 2 12 0
# 1 2 24 1
# 1 2 8 1
# 1 2 8 0
# 2 1 34 0
# 2 1 4 1
# 2 1 17 0
# 2 2 8 1
# 2 2 11 0
# 2 2 23 0
# 2 2 12 0
# 2 2 15 0
# 2 2 8 0
# 2 2 8 0
# 3 1 10 1
# 3 1 5 1 
# 3 2 25 0
# 3 2 8 1
# 3 2 11 0"
# melanoma = myDATstring.to.df(melan,m.names)
# write.csv(melanoma, file = "/Users/jamescutler/Desktop/Survival/melanoma.csv")
mel = read.csv("/Users/jamescutler/Desktop/Survival/melanoma.csv")

# And now the leukemia data:
leuk = read.csv("/Users/jamescutler/Desktop/Survival/leukemia.csv")
leuk.drug = leuk[leuk$group == 1,]
# setwd("/Users/jamescutler/Desktop/Survival/")
# write.csv(leuk.drug, file = "leukemia_drug.csv")

l.limits = c(6,12,24,36)
l.disc = 
  contToDisc(dataSet = leuk.drug,
             timeColumn = "time",
             intervalLimits = l.limits)
l.life = 
  lifeTable(l.disc,
            timeColumn = "timeDisc",
            censColumn = "delta",
            intervalBorders = myintBorders(l.limits))
l.life

myLTautoplot(l.life,l.limits)



# Myeloma:
myeloma = read.csv("/Users/jamescutler/Desktop/Survival/multiple_myeloma.csv")
m1 = myeloma[myeloma$protein == 1,]
m1
m0 = myeloma[myeloma$protein == 0,]
m0
m.limits = c(12,24,36,48,60,120)
m1.disc = contToDisc(dataSet = m1,
                     timeColumn = "time",
                     intervalLimits = m.limits)
m1.life = lifeTable(m1.disc,
                    timeColumn = "timeDisc",
                    censColumn = "status",
                    intervalBorders = myintBorders(m.limits))
m1.life

m0.disc = 
  contToDisc(m0,
             timeColumn = "time",
             intervalLimits = m.limits)
m0.life = lifeTable(m0.disc,
                    timeColumn = "timeDisc",
                    censColumn = "status",
                    intervalBorders = myintBorders(m.limits))
m0.life

m0.plot = myLTautoplot(m0.life,m.limits) 
m1.plot = myLTautoplot(m1.life,m.limits)

m0.plot +
  geom_line(data = m1.plot, mapping = aes(time,survival), col = "red") +
  geom_point(data = m1.plot, mapping = aes(time,survival), col = "red") +
  scale_colour_manual(name = "Protein",
                      values = c("red"="red","black"="black"),
                      labels = c("one","zero"))
myLTautoplot = function(disc.life,intBorders){
  new.disc.df = data.frame(time = intBorders,
                           survival = disc.life$Output$S)
  ggplot(new.disc.df, aes(time,survival)) +
    geom_line() +
    geom_point() +
    coord_cartesian(xlim = c(0,max(intBorders)),ylim = c(0,1)) +
    labs(title = "Life Table Survival Curve",
         x="Time",y="Survival Probability")
  #return(new.disc.df)
}
my.df = data.frame(protein_one = m1.life$Output$S,
                   protein_zero = m0.life$Output$S,
                   time = m.limits)
my.df
my.long = gather(my.df,
                 key = "protein",
                 value = "survival",
                 c("protein_one","protein_zero"))
ggplot(my.long, aes(x=time,y=survival,col=protein)) +
  geom_line() +
  geom_point() + 
  coord_cartesian(xlim = c(0,125),ylim = c(0,1)) +
  labs(title = "Life Table Survival Curve for Myeloma Data",
       x = "Time in months",y="Survival Probability")

# Now for the Nelson-Aalen cum hazard function plot:
m.fit = survfit(Surv(time,status) ~ protein, data = myeloma)
ggsurvplot(m.fit,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           fun = "cumhaz")

ggsurvplot(m.fit,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           fun = "event")


m.disc = 
  contToDisc(myeloma,
             timeColumn = "time",
             intervalLimits = m.limits)

m.life = 
  lifeTable(m.disc,
            timeColumn = "timeDisc",
            censColumn = "status",
            intervalBorders = myintBorders(m.limits))
m.life
m1.life
m0.life

bw = 35
reg0 = ksmooth(x=m.limits,y=m0.life$Output$hazard,kernel = "normal", bandwidth = bw)
plot(reg0$x,reg0$y,type = "l",lwd=2,col="blue",
     main = "Normal kernel-smoothed hazard functions",
     xlab = "Remission time (months)",
     ylab = "Estimated hazard rate")

reg1 = ksmooth(x=m.limits,y=m1.life$Output$hazard,kernel = "normal",bandwidth = bw)
points(reg1$x,reg1$y,type = "l",lwd=2,col="red")


myeloma.haz = data.frame(m0.haz = m0.life$Output$hazard,
                         m1.haz = m1.life$Output$hazard,
                         time = m.limits)
mhaz.long = 
  gather(myeloma.haz, 
         key = "protein",
         value = "hazard",
         c(m0.haz,m1.haz))

ggplot(mhaz.long, aes(x=time,y=hazard,col=protein)) +
  geom_smooth(method = "loess",level = .95) +
  scale_color_manual(values = c("blue","red")) +
  labs(title = "Loess-smoothed hazard functions",
       x="Remission time (months)",
       y="Estimated hazard rate")






######################################################################
######################################################################

# Week 4 lecture in-class data example

one = build_frame(
  "group", "y","delta"|
    "A"  , 4  , 1 |
    "A"  ,12  , 0 |
    "A"  ,15  , 1 |
    "A"  ,21  , 0 |
    "A"  ,23  , 1 |
    "B"  ,2   , 1 |
    "B"  ,6   , 0 |
    "B"  ,8   , 0 |
    "B"  ,10  , 1 |
    "B"  ,19  , 1
)

one.fit = survfit(Surv(y,delta) ~ group, data = one)
ggsurvplot(one.fit)

ggsurvplot(one.fit,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           fun = "event")
ggsurvplot(one.fit,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           fun = "cumhaz")
survdiff(Surv(y,delta) ~ group, data = one)
survdiff(Surv(y,delta) ~ group, data = one, rho = 0)


survdiff(Surv(time,status) ~ treatment, data = mel, rho = 0) # log-rank (Mantel-Haenszel)
survdiff(Surv(time,status) ~ treatment, data = mel, rho = 1) # Peto & Peto mod. of Gehan-Wilcoxon


mel.fit = survfit(Surv(time,status) ~ treatment, data = mel)
mel.age.fit = survfit(Surv(time,status) ~ age, data = mel)
ggsurvplot(mel.age.fit)

m.dt = data.table(mel)
m.dt[age == 2,][order(time)]

mel1 = mel[mel$treatment == 1,]
mel1.age.fit = survfit(Surv(time,status) ~ age, data = mel1)
ggsurvplot(mel1.age.fit)











###########################################################################

# HW Week 4

1-pchisq(3.515,df=1)
qnorm(.975)
qchisq(.95,df=1)

-4.565/2.435



canc = "1 23 1
1 47 1
1 69 1
1 70 0
1 71 0
1 100 0
1 101 0
1 148 1
1 181 1
1 198 0
1 208 0
1 212 0
1 224 0
2 5 1
2 8 1
2 10 1
2 13 1
2 18 1
2 24 1
2 26 1
2 26 1
2 31 1
2 35 1
2 40 1
2 41 1
2 48 1
2 50 1
2 59 1
2 61 1
2 68 1
2 71 1
2 76 0
2 105 0
2 107 0
2 109 0
2 113 1
2 116 0
2 118 1
2 143 1
2 154 0
2 162 0
2 188 0
2 212 0
2 217 0
2 225 0"
canc = myDATstring.to.df(canc,"stain time status")
head(canc)

cdt = data.table(canc)

# Step 1: order it be time
ansnew=cdt[order(time)]

# Step 2: create a new df with just events and no censored
ans = ansnew[status == 1]

# Step 3: make sure there are no duplicate times, by summing events by time 
## (leaving you with the right number of rows for the data frame)
ans2 = ans[,sum(status),by="time"]
colnames(ans2) = c("time","deaths")
# Getting the matching right is key, and it makes sense:
ans2$stain = ans$stain[match(ans2$time,ans$time)]
ans2 = ans2[,c(2,3,1)]
ans2

# Step 4: add a d1j column
ans2$d1j = NA
deaths1 = ans2[stain == 1,which=TRUE]
ans2$d1j[deaths1] = ans2$deaths[deaths1]
ans2$d1j[-deaths1] = 0
ans2

# Step 5: add an n1j column (which requires adding a d1j to the full data frame)
## Add a column indicating stain-1-only censored people, to the full data frame:
cens1 = ansnew[stain == 1 & status == 0,which=TRUE]
ansnew$n1j.cens = NA
ansnew$n1j.cens[cens1] = 1
ansnew$n1j.cens[-cens1] = 0
ansnew
## Add a d1j column to the full data frame:
ansnew$d1j = NA
deaths1 = ansnew[stain == 1 & status == 1,which=TRUE]
ansnew$d1j[deaths1] = 1
ansnew$d1j[-deaths1] = 0
ansnew
## The genius code for creating n1j, first for the full data frame:
n1j = vector()
for (i in 1:nrow(ansnew)){
  n1j[i] = ansnew[stain == 1,.N] - 
    ( sum(ansnew$d1j[1:i-1]) + sum(ansnew$n1j.cens[1:(i-1)]) )
}
n1j
ansnew$n1j = n1j
ansnew
## And now match an n1j column to the log-rank table:
ans2$n1j = 
  ansnew$n1j[match(ans2$time,ansnew$time)]
ans2

# Step 6 (REPEAT STEP 4): add a d2j column
ans2$d2j = NA
deaths2 = ans2[stain == 2,which=TRUE]
ans2$d2j[deaths2] = ans2$deaths[deaths2]
ans2$d2j[-deaths2] = 0
ans2

# Step 7 (REPEAT STEP 5): add an n2j column
## Add a stain-2-only censored column to the full data frame:
cens2 = ansnew[stain == 2 & status == 0,which=TRUE]
ansnew$n2j.cens = NA
ansnew$n2j.cens[cens2] = 1
ansnew$n2j.cens[-cens2] = 0
ansnew
## Add a d2j column to the full data frame:
ansnew$d2j = NA
fulldeaths2 = ansnew[stain == 2 & status == 1,which=TRUE]
ansnew$d2j[fulldeaths2] = 1
ansnew$d2j[-fulldeaths2] = 0
ansnew
## The genius code for creating n2j, first for the full data frame:
n2j = vector()
for (i in 1:nrow(ansnew)){
  n2j[i] = ansnew[stain == 2,.N] - 
    ( sum(ansnew$d2j[1:i-1]) + sum(ansnew$n2j.cens[1:(i-1)]) )
}
n2j
ansnew$n2j = n2j
ansnew
## And now match an n1j column to the log-rank table:
ans2$n2j = 
  ansnew$n2j[match(ans2$time,ansnew$time)]
ans2

# Step 8: add d1j and d2j to get dj
ans2$dj = ans2$d1 + ans2$d2j         # I have no idea why this works, but it does, even though sum(ans2$d1j,ans2$d2j) doesn't
ans2

# Step 9: add n1j and n2j to get nj
ans2$nj = ans2$n1j + ans2$n2j
ans2

# Step 10: get e1j
ans2$e1j = ans2$n1j*ans2$dj/ans2$nj
ans2

# Step 11: get v1j
ans2$v1j = ( ans2$n1j*ans2$n2j*ans2$dj*(ans2$nj-ans2$dj) ) /
  ( (ans2$nj^2)*(ans2$nj - 1) )
ans2

# Step 12: get O1j - Ej column
ans2$OminusE = ans2$d1j - ans2$e1j
ans2

# Now sum OminusE, square the sum, and divide by sum(v1j) to get your X^2 statistic:
X2L = (sum(ans2$OminusE)^2)/sum(ans2$v1j)
X2L.pvalue =
  1-pchisq(X2L,df=1)

# Or, sum OminusE, and divide by sqrt(sum(v1j)) to get the XL statistic:
XL = sum(ans2$OminusE)/sqrt(sum(ans2$v1j))
abs(XL)
qnorm(.975)









##########################################################################
##########################################################################

# Week 5 - Cox models

# Libraries required for Cox models in R:
## survival
## survminer

myel = read.csv("/Users/jamescutler/Desktop/Survival/multiple_myeloma.csv")
head(myel)


# For running a bunch of UNIVARIATE analyses:
covariates <- c("age","sex","bun","ca","hb","pcells","protein")
univ_formulas <- sapply(covariates,
                        function(x) {as.formula(paste('Surv(time, status)~', x))}
)
univ_models <- lapply(univ_formulas, function(x) {coxph(x, data = myel)}
)
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)


# For a multivariate Cox regression analysis:
res.cox = coxph(Surv(time,status) ~ age + sex + bun + ca + hb + pcells + protein,
                data = myel)
summary(res.cox)


# Test for bun*hb interaction:
res.cox = coxph(Surv(time,status) ~ 
                  bun + hb + bun*hb,
                data = myel)
summary(res.cox)




##########

# bladder cancer DAT to csv:
# blad = [bunch of data]
# setwd("/Users/jamescutler/Desktop/Survival/")
# write.csv(blad,file = "bladder_cancer.csv")





##########################################################################
##########################################################################

# Week 5 HW

# Datasets
## bladder cancer
## prostatic cancer
pros = read.csv("/Users/jamescutler/Desktop/Survival/Prostate_Cancer.csv")
blad = read.csv("/Users/jamescutler/Desktop/Survival/bladder_cancer.csv")

head(blad)
# blad$treat = as.factor(blad$treat)
# blad$status = as.factor(blad$status)
# summary(blad)



blad$init.group = NA
blad$init.group[blad$init == 1] = 1
blad$init.group[blad$init %in% 2:3] = 2
blad$init.group[which(is.na(blad$init.group))] = 3
blad$init.group
hist(blad$init.group)

bdt = data.table(blad)
hist(blad$init)
bdt[,.N,by="init"]

# Test to see whether the three hazard functions for init.group=1,2,3 are different.

placebo = bdt[treat == 1]
drug = bdt[treat == 2]

# Log-rank test on the "survival curves" (same as doing it on hazard functions)
s.diff = survdiff(Surv(time,status) ~ init.group, data = placebo)
1-pchisq(s.diff$chisq,df=2)     # The p-value with more decimal places

drug.diff = survdiff(Surv(time,status) ~ init.group, data = drug)
1-pchisq(drug.diff$chisq,df=2)  # The p-value with more decimal places
drug.fit = survfit(Surv(time,status) ~ init.group, data = drug)
ggsurvplot(drug.fit,
           conf.int = TRUE) +
  labs(title = "Survival curves for the three groups of initial number of tumors, in the drug arm",
       x="Time in months")
ggsurvplot(drug.fit,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           fun = "event")

# Log-rank test for trend using package survMisc's comp() function
p.fit = survfit(Surv(time,status) ~ init.group, data = placebo)
summary(p.fit)
# p.trend.test = comp(p.fit)
p.ten = ten(p.fit)
comp(p.ten)
1-pnorm(1.08053)

# Now trend for drug arm:
drug.ten = ten(drug.fit)
comp(drug.ten)
1-pnorm(3.0764)
# autoplot(p.fit) - Bizarrely, you can't use it anymore once you call survMisc
ggsurvplot(p.fit, 
           conf.int = TRUE) +
  labs(title = "Survival curves for the three groups of initial number of tumors, in the placebo arm",
       x="Time in months")


# muhaz for hazard plots
p.haz = muhaz(placebo$time,placebo$status)
plot(p.haz)
library(splines)
x<-1:5
y <- c(0.31, 0.45, 0.84, 0.43, 0.25)
yy <-predict(interpSpline(x, y))
plot(x, y)
lines(yy)


blad.cox = coxph(Surv(time,status) ~ init + size,
                 data = drug)
summary(blad.cox)

blad.coxplace = coxph(Surv(time,status) ~ init + size,
                      data = placebo)
summary(blad.coxplace)




# Part 2: prostate cancer and Cox regression models
pros.cox = coxph(Surv(time,status) ~ treatment + age + shb + size + index, 
                 data = pros)
summary(pros.cox)

tidy(pros.cox)

# Awesome:
ggforest(pros.cox, data = pros)


# Part 2 b:
# pros$tx = NA
# pros$tx[pros$treatment == 1] = 0
# pros$tx[pros$treatment == 2] = 1
# Should just do:
pros$tx = pros$treatment - 1

# Main effects with just tx, size, and index
pros.coxb = coxph(Surv(time,status) ~ treatment + size + index,
                  data = pros)
summary(pros.coxb)

# Full model - with tx:size and tx:index interactions added
pros.coxb2 = coxph(Surv(time,status) ~ treatment + size + index +
                     treatment:size + treatment:index, 
                   data = pros)
summary(pros.coxb2)

tidy(pros.coxb2)
ggforest(pros.coxb2, data = pros)


# The following uses the SAS-derived X^2_Cox test statistic and DF
1-pchisq(2.867,2)                # P=0.23847

# And this is the result from R
anova(pros.coxb,pros.coxb2)      # P=0.2384
# THEY'RE THE SAME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Conclusion: We fail to reject the null and conclude that the reduced model is adequate.



##########################################################################
##########################################################################

# Week 6 HW

# anderson = [bunch of data]
# setwd("/Users/jamescutler/Desktop/Survival/")
# write.csv(leuk, file = "Anderson_leukemia.csv")
leuk = read.csv("/Users/jamescutler/Desktop/Survival/Anderson_leukemia.csv")

# This is the right cox model, and thus will be used in the creation of the four-curve plot
l.cox = coxph(Surv(time,delta) ~ sex + logWBC + trt +
                trt:sex,
              data = leuk)
summary(l.cox)

# This one removes all the coefficients except for the one for logWBC!
# and.cox = coxph(Surv(time,delta) ~ strata(sex) + logWBC + strata(trt),
#                 data = leuk)
# summary(and.cox)

and.cox = coxph(Surv(time,delta) ~ sex + logWBC + trt + trt:sex)

# Unfortunately, I have no idea how to get the hazard ratios for 
## (drug, female) and (drug, male).
# Here's an attempt though:
ldt = data.table(leuk)
leuk$trt.male = NA
leuk$trt.female = NA
leuk$trt.male[ldt[sex == 0 & trt == 0,which=TRUE]] = "M,P"
leuk$trt.male[ldt[sex == 0 & trt == 1,which=TRUE]] = "M,tD"
leuk$trt.female[ldt[sex == 1 & trt == 0,which=TRUE]] = "F,P"
leuk$trt.female[ldt[sex == 1 & trt == 1,which=TRUE]] = "F,tD"
leuk.male = leuk[-which(is.na(leuk$trt.male)),]
leuk.female = leuk[-which(is.na(leuk$trt.female)),]
# Attempt number 1 (not conditioning on logWBC):
l.cox.male = coxph(Surv(time,delta) ~ trt.male,
                   data = leuk.male)
l.cox.female = coxph(Surv(time,delta) ~ trt.female,
                     data = leuk.female)
summary(l.cox.male)
summary(l.cox.female)
## They're fairly close to the SAS results.
# Attempt number 2 (holding logWBC fixed):
l.cox.male = coxph(Surv(time,delta) ~ logWBC + trt.male,
                   data = leuk.male)
l.cox.female = coxph(Surv(time,delta) ~ logWBC + trt.female,
                     data = leuk.female)
summary(l.cox.male)
summary(l.cox.female)
## They're still sort of close, but it actually gets less like the SAS results
## after conditioning on logWBC. Weird.


# A look at how the strata() function works (even though I don't need it here):
a <- factor(rep(1:3,4), labels=c("low", "medium", "high"))
b <- factor(rep(1:4,3))
levels(strata(b))
levels(strata(a,b,shortlabel=TRUE))
matrix(strata(leuk$sex,leuk$trt))

# Get four curves by creating a dataframe with four new rows:
sex.trt = with(leuk,
               data.frame(trt = c(1,1,0,0),
                          sex = c(1,0,1,0),
                          logWBC = rep(2,4) 
               )
)
sex.trt

# This is the right one!!!:
## Source: http://www.sthda.com/english/wiki/cox-proportional-hazards-model
## EXTREMELY IMPORTANT INFORMATION!!!!!!!!!! ^
leuk.fit = survfit(l.cox, data = leuk, newdata = sex.trt)
ggsurvplot(leuk.fit,
           conf.int = FALSE,
           legend.labs = c("Drug, Female",
                           "Drug, Male",
                           "Placebo, Female",
                           "Placebo, Male")) +
  labs(title = "Survival Functions by Treatment Group and Gender",
       x="Survival time (weeks)",y="Survival Function")



# Part 2 - melanoma
mel = read.csv("/Users/jamescutler/Desktop/Survival/melanoma.csv")
mel

# Age as a categorical variable:
mel$age.groups = NA
mel$age.groups[mel$age == 1] = "21-40"
mel$age.groups[mel$age == 2] = "41-60"
mel$age.groups[mel$age == 3] = "61-"
mel$age.groups

m.cox123 = coxph(Surv(time,status) ~ age.groups,
                 data = mel)
summary(m.cox123)
# You get hazard ratios for 41-60 and 61+ that are both relative to the 21-40 baseline.
## This approach reveals a depression in hazard in the middel age group, and an increase
## in the oldest age group (non-linear relationship between age and hazard?).

# Age as a continuous variable:
m.cox123.cont = coxph(Surv(time,status) ~ age,
                      data = mel)
summary(m.cox123.cont)
# You get one hazard ratio telling you how much hazard changes with one unit increase in 
## age (1 to 2, 2 to 3). 
## This approach conceals the complexity revealed in the categorical approach.


# Compare the categorical model to the continuous model
## Our question is, Is treating age as a continuous covariate adequate?
anova(m.cox123,m.cox123.cont) # The p-value is equivalent to: 1-pchisq(.9092,1)
anova(m.cox123.cont,m.cox123) # Same result (order of arguments doesn't matter)
## However, the survival package documentation does say that "It is conventional 
## to list the models from smallest to largest, but this is up to the user."
## Note: The test statistics and p values compare the reduction in log-likelihood 
## for each row. 

(29.669-29.215)*2

1-pchisq(.317,1)

1-pchisq(3.237,5)







##################################################################################

# pbc = [bunch of data]
# cirrhosis = 
#   myDATstring.to.df(pbc,"id futime status drug age sex ascites hepato spiders edema bili chol albumin copper alk_phos sgot trig platelet protime stage")
# mynames = 
# "stage id futime status drug age sex ascites hepato spiders edema bili chol albumin copper alk_phos sgot trig platelet protime"             
# mynames = unlist(strsplit(mynames," "))
# 
# one = cirrhosis[1,2]
# almost = cirrhosis[2:(nrow(cirrhosis)-1),]
# head(almost)
# colnames(almost) = mynames
# head(almost)
# almost = almost[,c(2:20,1)]
# tail(almost)
# head(almost)
# almost$stage = 
# c(almost$stage[2:nrow(almost)],4)
# head(almost)
# tail(almost)
# 
# pbc = rbind(one,
#             almost)
# head(pbc)
# tail(pbc)
# 
# # Now I can finally write it as a csv
# write.csv(pbc, file = "/Users/jamescutler/Desktop/Survival/pbcirrhosis.csv")
pbc = read.csv("/Users/jamescutler/Desktop/Survival/pbcirrhosis.csv")

pbc = pbc[1:312,]
pbc$status
pbc$status[which(pbc$status == 2)] = 1
pbc$status
pbc$trt = 2 - pbc$drug
pbc$trt
pbc$drug
pbc$yrs.age = pbc$age/365
pbc$yrs.age

p.cox = coxph(Surv(futime,status) ~ 
                sex + 
                yrs.age + 
                edema + 
                trt +
                bili +
                strata(stage),
              data = pbc)

ggsurvplot(survfit(p.cox, data = pbc),
           fun = "cumhaz")





# In class exercise on addicts.dat dataset
# addict = [bunch of data]
# write.csv(ad, file = "/Users/jamescutler/Desktop/Survival/addicts.csv")
addict = read.csv("/Users/jamescutler/Desktop/Survival/addicts.csv")

a.cox = coxph(Surv(time,status) ~ 
                clinic + 
                prison +
                dose,
              data = ad)
ggsurvplot(survfit(a.cox,data=ad),
           fun = "cumhaz")

a.cox2 = coxph(Surv(time,status) ~ 
                 strata(clinic) +
                 prison +
                 dose,
               data = ad)
ggsurvplot(survfit(a.cox2,data=ad),
           fun = "cumhaz")
summary(a.cox)
summary(a.cox2)
extractAIC(a.cox)      # My AIC value is almost exactly the same as SAS. (1352)
extractAIC(a.cox2)     # My AIC value is almost exactly the same as SAS. Heck yeah! (1199)


# Proportional hazards assumption is simply the assumption that the hazards 
## are proportional (log of curves are parallel). 
# If the assumption is rejected, then we reject [dummy variable stratification]
## and we do true stratification.
cox.zph(a.cox)     # This gets me a p value that is significant like SAS, but lower (?).
cox.zph(a.cox2)    # This is probably not the right model to run the test on.







#################################################################################
#################################################################################

# Week 7 HW

addict = read.csv("/Users/jamescutler/Desktop/Survival/addicts.csv")

# Making prison categorical changes nothing
# addict$prison[addict$prison == 0] = "no"
# addict$prison[addict$prison == 1] = "yes"
# addict$prison = as.factor(addict$prison)

# The easy way--typing factor(prison) inside the model--does nothing either

# Stratify clinic (true stratification) - use exact method for ties
a.cox = coxph(Surv(time,status) ~ factor(prison) + dose +
                strata(clinic),
              data = addict, ties = "exact")
# Plot it
a.fit = survfit(a.cox,data=addict)
ggsurvplot(a.fit,
           linetype = "strata") +
  labs(title = "Cox model of survival probability with true stratification of clinic",
       x="Time in days")

# Don't stratify clinic (for the sake of contrast with above) - use exact method for ties
a.cox2 = coxph(Surv(time,status) ~ factor(prison) + dose +
                 + clinic,
               data = addict, ties = "exact")
a.fit2 = survfit(a.cox2,data=addict)
ggsurvplot(a.fit2) +
  labs(title = "Cox model of survival probability without stratification of clinic")

# What is the effect of prison and dose when clinic is true stratified
summary(a.cox)


# Use efron method for ties to compare to the results for prison and dose using the exact method
a.coxEfron = coxph(Surv(time,status) ~ factor(prison) + dose +
                     strata(clinic),
                   data = addict, method = "efron")
summary(a.coxEfron)
summary(a.cox)


# Part 2

leu = read.csv("/Users/jamescutler/Desktop/Survival/leukemia.csv")

# Change 0 to control and 1 to drug in variable "group"
leu$group[leu$group == 0] = "Control"
leu$group[leu$group == 1] = "Drug"

l.nonstrata = coxph(Surv(time,delta) ~ group,
                    data = leu)
summary(l.nonstrata)
l.nonstrata.fit = survfit(l.nonstrata)
ggsurvplot(l.nonstrata.fit, data = leu)

# Stratified Cox model and ordinary KM survift object
l.cox = coxph(Surv(time,delta) ~ strata(group),
              data = leu)
l.km = survfit(Surv(time,delta) ~ group, 
               data = leu)

l.coxfit = survfit(l.cox)

# Combine the two fits into a list
comb.fit = list(Cox = l.coxfit, KM = l.km)

# Plug the list object into ggsurvplot_combine
## Source: https://rpkgs.datanovia.com/survminer/reference/ggsurvplot_combine.html
## For linetype = "strata", see http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
ggsurvplot_combine(comb.fit,data=leu,
                   linetype = "strata") +
  labs(title = "KM survival curves and Cox model estimates of survival curves for drug and control groups",
       x="Remission time in weeks")
# Voila

# Dummy variable cox model vs KM curves
## There is no way to do a cox model on a dataset where there is nothing but one
## one-level covariate. The model will be empty. So what the heck is dummy variable
## stratification?



and.cox = coxph(Surv(time,delta) ~ 
                  logWBC
                +
                  trt*sex
                ,
                data = leuk)
summary(and.cox)
#################################################################################
#################################################################################

# Week 8 in class

pbc = read.csv("/Users/jamescutler/Desktop/Survival/pbcirrhosis.csv")
pbc = pbc[1:312,]
pbc$status[pbc$status == 1] = 0
pbc$status[pbc$status == 2] = 1



head(pbc)
pbc$treat = 2 - pbc$drug
pbc$age.year = pbc$age/365
pbc$logalb = log(pbc$albumin)
pbc$logbili = log(pbc$bili)



b.cox = coxph(Surv(futime,status) ~ bili, data = pbc)
plot(pbc$bili,b.cox$residuals,ylim = c(-1,1))
range(pbc$bili)


# Untransformed
ggcoxfunctional(Surv(futime,status) ~ bili + albumin + age.year + edema, 
                data = pbc,
                ylim = c(-1,1))

# Log transformed - much more linear!
pbc.cox = coxph(Surv(futime,status) ~ log(bili) + log(albumin) + age.year + edema,
                data = pbc)
ggcoxfunctional(pbc.cox, 
                data = pbc,
                ylim = c(-2,1))

# Test proportional hazards assumption
test.ph = cox.zph(pbc.cox)
test.ph

# Plot martingale residuals vs linear predictors
plot(pbc.cox$linear.predictors,pbc.cox$residuals, ylim = c(-3,1))
abline(h=0,col = "gray")

# Plot deviance residuals vs linear predictors
deviance.res = 
  resid(pbc.cox, type = "deviance")
plot(pbc.cox$linear.predictors,deviance.res, ylim = c(-3,3))
abline(h=0,col = "gray")

# Plot deviance residuals (full model) by each (transformed) covariate










#####################################################################################
#####################################################################################

# Week 10 - checking the Cox model, continued

ggcoxzph(test.ph)

# dfbeta vs dfbetas (standardized)
ggcoxdiagnostics(pbc.cox, type = "dfbeta",
                 linear.predictions = FALSE)

ggcoxdiagnostics(pbc.cox, type = "dfbetas",
                 linear.predictions = FALSE)

# deviance
ggcoxdiagnostics(pbc.cox, type = "deviance",
                 linear.predictions = FALSE)


# Taking continuous data and making it categorical using cut()
cut(pbc$bili, c(0,1.45,3.25,6.75,max(pbc$bili)))                   # AMAZING!!!
bili.strat =
  cut(pbc$bili, c(0,1.45,3.25,6.75,max(pbc$bili)+.1), right = FALSE) # AMAZING!!! - you can have it include the lower, exlcude the upper, or vice versa
summary(bili.strat)                                                # AMAZING!!!
barplot(summary(bili.strat))                                       # AMAZING!!! yeah!

alb.strat = cut(pbc$albumin, c(0,3.1,3.4,3.67,max(pbc$albumin)+.1), right = FALSE)
summary(alb.strat)
barplot(summary(alb.strat))

age.strat = cut(pbc$age.year, c(0,46,53,61,max(pbc$age.year)+.1), right = FALSE)
summary(age.strat)
barplot(summary(age.strat))

pbc$biliStrat = bili.strat
pbc$albStrat = alb.strat
pbc$ageStrat = age.strat




# In class exercise
leu = read.csv("/Users/jamescutler/Desktop/Survival/Anderson_leukemia.csv")
head(leu)
leu$wbc = exp(leu$logWBC)
head(leu)

l.wbc = coxph(Surv(time,delta) ~ wbc, data = leu)

# Martingale residuals by predictors
plot(l.wbc$linear.predictors,l.wbc$residuals)
abline(h=0,col="gray")
plot(leu$wbc,l.wbc$residuals)                 # Same exact positions as above, but with different x axis

# Formal test for PH assumption
test.wbc = cox.zph(l.wbc)
test.wbc
plot(test.wbc)
abline(h=0,col="red",lty = 2)

# Formal test for functional form
ggcoxfunctional(l.wbc, data=leu)              # Huge non-linearity
# install.packages("goftte")
library(goftte)
fcov(l.wbc)
.135/2
fcov(l.wbc,type.test = "Liu")
# BIG TAKEAWAY HERE: If you have a formal test result (P-value barely not less than .05) 
## that suggests that there is no problem with the linearity, but your dataset is small,
## then you might have a problem with the test being underpowered--it might not be powerful
## enough to detect the non-linearity that is in fact there! That's why in this case the 
## graphical tool is more trustworthy than the 0.135, or 0.061 [0.073 in SAS] p-value!

# Now compare to logWBC
log.wbc = coxph(Surv(time,delta) ~ logWBC, data = leu)
# Formal test for PH assumption
test.logwbc = cox.zph(log.wbc)
test.logwbc
# Formal test for functional form
ggcoxfunctional(log.wbc, data = leu)    # Not nearly as much non-linearity
# Martingale resituals by predictors
plot(log.wbc$linear.predictors,log.wbc$residuals); abline(h=0,col="red")
# Martingale by original data
plot(leu$wbc,log.wbc$residuals); abline(h=0,col="red")


# Test PH assumption for Trt and logWBC
meMod = coxph(Surv(time,delta) ~ trt + logWBC, data = leu)
test.me = cox.zph(meMod)
test.me

# Test functional form for Trt and logWBC
ggcoxfunctional(meMod, data = leu)

# Deviance residuals by linear predictor - SHOWS JUST ONE PLOT
ggcoxdiagnostics(meMod, type = "deviance",
                 linear.predictions = TRUE)  # could it be one plot because of this? No. It's always one plot. 

# Index plot of standardized delta betas - SHOWS TRT AND LOGWBC
ggcoxdiagnostics(meMod, type = "dfbetas",
                 linear.predictions = FALSE) # could it be two plots because of this? No. It's supposed to be false. 
## The warning says linear.predictions only works for martingales and deviance

ggcoxdiagnostics(meMod, type = "dfbeta",
                 linear.predictions = FALSE)













# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
################### LOAD mm DATA FROM WEEK 10 HW BELOW ###################
# ALL 21 MODELS WITH A TWO-WAY INTERACTION
covariates <- c("age","sex","bun","ca","hb","pcells","protein")

# Main effects
part1 = paste(covariates,collapse = "+")

# All possible combinations of 2
all.ints = combn(covariates,2)

# Paste the combinations together into interactions with ':'
part2 = character()
for (i in 1:(length(all.ints)/2)){
  part2[i] = paste(all.ints[(i*2 - 1):(i*2)],collapse = ":")
}
part2

# Create the full list of terms to go in the formula
theterms = paste(part1,part2,sep = "+")

# Create the formulas - COULDN'T IT JUST BE int.formulas = as.formulas(paste('Surv(time,status)~',theterms)) ???
int.formulas =
  sapply(theterms, function(x){
    as.formula(paste('Surv(time,status)~',x))
  })

# Create the cox models
int.models =
  lapply(int.formulas, function(x){
    coxph(x,data=mm)
  })
summary(int.models[[1]])

# Print the summaries of the models
for (i in 1:length(int.models)){
  print(summary(int.models[[i]]))
}

# The only interactions that are significant are:
## age:sex 
## pcells:protein
summary(int.models[[1]])
summary(int.models[[21]])
#####################################################################################
#####################################################################################

# Week 10 HW
mm = read.csv("/Users/jamescutler/Desktop/Survival/multiple_myeloma.csv")
head(mm)
mm = mm[,2:11]
mmdt = data.table(mm)

ggsurvplot(survfit(Surv(time,status) ~ sex, data = mm),
           surv.median.line = "hv"
)
# Univariate analyses
covariates <- c("age","sex","bun","ca","hb","pcells","protein")

univ_formulas <- sapply(covariates, function(x){
  as.formula(paste('Surv(time, status)~', x))
})

univ_models <- lapply(univ_formulas, function(x){
  coxph(x, data = mm)
})

for (i in 1:length(univ_models)){
  print(summary(univ_models[[i]]))
}

# The only IVs significantly associated with survival time in the univariate
## analyses are:
## bun 
## hb
summary(univ_models[[3]])
summary(univ_models[[5]])


# Bence-Jones protein associated with longer median survival (weird)
ggsurvplot(survfit(Surv(time,status) ~ protein, data = mm),
           surv.median.line = "hv",
           break.time.by = 10) 


# Descriptive statistics 
## SEE Survival_HW10_descriptive.Rmd


# Functional form formal tests - FOR CONTINUOUS ONLY
fcov(univ_models[[1]], type.test = "Liu")
# fcov(univ_models[[2]], type.test = "Liu")       # <.001 - not continuous though
fcov(univ_models[[3]], type.test = "Liu")
fcov(univ_models[[4]], type.test = "Liu")
fcov(univ_models[[5]], type.test = "Liu")
fcov(univ_models[[6]], type.test = "Liu")
# fcov(univ_models[[7]], type.test = "Liu")       # < .001 - not continuous though

# Functional form plots 
## THESE WERE THE SAME AS SAS OUTPUT IN THE WEEK 10 IN CLASS EXERCISE, BUT
## NOW THEY'RE DIFFERENT (???)
ggcoxfunctional(univ_models[[1]], data = mm)
# ggcoxfunctional(univ_models[[2]], data = mm)      # NOT CONTINUOUS! DON'T NEED TO ASSESS FUNCTIONAL FORM
ggcoxfunctional(univ_models[[3]], data = mm)
ggcoxfunctional(univ_models[[4]], data = mm)
ggcoxfunctional(univ_models[[5]], data = mm)
ggcoxfunctional(univ_models[[6]], data = mm)
# ggcoxfunctional(univ_models[[7]], data = mm)  # NOT CONTINUOUS! DON'T NEED TO ASSESS FUNCTIONAL FORM

# Martingale residual plots
par(mfrow=c(2,3))
plot(mm$age,resid(univ_models[[1]]))
lines(lowess(resid(univ_models[[1]]) ~ mm$age, f=2/3), col = "red")
plot(mm$bun,resid(univ_models[[3]]))
lines(lowess(resid(univ_models[[3]]) ~ mm$bun, f=2/3), col = "red")
plot(mm$ca,resid(univ_models[[4]]))
lines(lowess(resid(univ_models[[4]]) ~ mm$ca, f=2/3), col = "red")
plot(mm$hb,resid(univ_models[[5]]))
lines(lowess(resid(univ_models[[5]]) ~ mm$hb, f=2/3), col = "red")
plot(mm$pcells,resid(univ_models[[6]]))
lines(lowess(resid(univ_models[[6]]) ~ mm$pcells, f=2/3), col = "red")
par(mfrow=c(1,1))

# ... of just bun and hb
par(mfrow=c(1,2))
plot(mm$bun,resid(univ_models[[3]]),
     main = "Martingale residuals for blood urea nitrogen",
     xlab = "Blood urea nitrogen (mg/dL)",
     ylab = "Martingale residuals")
lines(lowess(resid(univ_models[[3]]) ~ mm$bun, f=2/3), col = "red")
plot(mm$hb,resid(univ_models[[5]]),
     main = "Martingale residuals for hemoglobin",
     xlab = "Hemoglobin (g/dL)",
     ylab = "Martingale residuals")
lines(lowess(resid(univ_models[[5]]) ~ mm$hb, f=2/3), col = "red")
par(mfrow=c(1,1))

# No need to do transformations, right?

# Main effects multivariable model - WITH ALL 7
mod.ma = coxph(Surv(time,status) ~ age+sex+bun+ca+hb+pcells+protein, data = mm)
summary(mod.ma)
# Only bun is significant. hb is marginally significant.
## But SAS has P=0.049 for hb, instead of P=0.052

# Main effects multivariable mode - WITH JUST BUN AND HB
mm$`Blood urea nitrogen (mg/dL)` = mm$bun
mm$`Hemoglobin (g/dL)` = mm$hb
mod.bunhb = coxph(Surv(time,status) ~ 
                    `Blood urea nitrogen (mg/dL)` + 
                    `Hemoglobin (g/dL)`, data = mm)
summary(mod.bunhb)
# Both are significant - THIS IS THE MODEL TO USE

# Interaction model - bun + hb + bun:hb
mod.bunhb.int = coxph(Surv(time,status) ~ bun*hb, data = mm)
summary(mod.bunhb.int)
# Nothing is significant




# Now do LRT comparing interaction model to just main effects 
anova(mod.bunhb,mod.bunhb.int)
# The reduced model should be accepted.




# Functional form for covariates in mod.bunhb
ggcoxfunctional(Surv(time,status) ~ 
                  `Blood urea nitrogen (mg/dL)` + `Hemoglobin (g/dL)`,  
                data = mm,
                title = "Martingal residuals from multivariable cox model")
fcov(mod.bunhb, type.test = "Liu")


# PH assumption
test.bunhb = cox.zph(mod.bunhb)
test.bunhb
ggcoxzph(test.bunhb)


# Outliers and influential observations
ggcoxdiagnostics(mod.bunhb, type = "dfbeta",
                 linear.predictions = FALSE)

ggcoxdiagnostics(mod.bunhb, type = "deviance",
                 linear.predictions = FALSE)

hist(mm$bun)
bun.m = mean(mm$bun)
bun.sd = sd(mm$bun)
bun.m + 2*bun.sd
max(mm$bun)
sort(mm$bun)
hist(mm$hb)


# Kaplan Meier or cox model survival or hazard curves?
mod.sex.strat = coxph(Surv(time,status) ~ strata(sex), data = mm)
summary(mod.sex.strat) # FASCINATING! YOU CAN'T JUST STRATIFY BY ONE VARIABLE!
splot.sex.strat = survfit(mod.sex.strat, data = mm)
ggsurvplot(splot.sex.strat,
           linetype = "strata")

mod.protein.strat = coxph(Surv(time,status) ~ strata(protein), data = mm)
ggsurvplot(survfit(mod.protein.strat, data = mm),
           linetype = "strata")





# Doing the dfbeta plots myself
mydfbeta = residuals(mod.bunhb, type = "dfbeta")
beta.bun = mydfbeta[,1]
beta.hb = mydfbeta[,2]
plot(1:48,beta.bun)
plot(1:48,beta.hb)
plot(mm$bun,beta.bun)
abline(h=c(-0.005816,0.005816), col = "grey")
plot(mm$hb,beta.hb)
# Get the SE of the beta coefficients for bun and for hb from the model, and 
## check that none of the deltas exceed it in their distance from zero. 








################################################################################
################################################################################

# Week 12

# Needs library(flexsurv)
# install.packages("flexsurv")
# library(flexsurv)
leu = read.csv("/Users/jamescutler/Desktop/Survival/leukemia.csv")
head(leu)
leu = leu[,2:4]
head(leu)


# Create the models
l.exp = survreg(Surv(time,delta) ~ group, data = leu, dist = "exponential")
summary(l.exp)
l.weib = survreg(Surv(time,delta) ~ group, data = leu, dist = "weibull")
summary(l.weib)
l.lognorm = survreg(Surv(time,delta) ~ group, data = leu, dist = "lognormal")
summary(l.lognorm)
# Note, this one requires flexsurv::flexsurvreg
l.gamma = flexsurvreg(Surv(time,delta) ~ group, data = leu, dist = "gamma")
l.gamma



# Now make the comparisons via LRT
myLRT = function(loglik1,loglik2, df1,df2){
  neg2logL1 = -2*loglik1
  neg2logL2 = -2*loglik2
  lrt = abs(neg2logL1-neg2logL2)
  finaldf = abs(df1-df2)
  if (finaldf == 0){
    finaldf = 1
    pval = 1-pchisq(lrt,finaldf)
  }else{
    pval = 1-pchisq(lrt,finaldf)
  }
  print(neg2logL1)
  print(neg2logL2)
  print("LRT ... DF ... P-value")
  return(c(lrt,finaldf,pval))
}
myLRT(l.exp$loglik[2],l.weib$loglik[2], 1,2); anova(l.exp,l.weib) # accept Weibull?
myLRT(l.exp$loglik[2],l.gamma$loglik, 1,3)       # P=0.124  vs .12; accept exponential
myLRT(l.weib$loglik[2],l.gamma$loglik, 2,3)      # P=0.5996 vs .54; accept Weibull
myLRT(l.lognorm$loglik[2],l.gamma$loglik, 2,3)   # P=0.468  vs .43; accept lognormal



# Now fit a log-logistic
l.loglogis = survreg(Surv(time,delta) ~ group, data = leu, dist = "loglogistic")
summary(l.loglogis)



# And compare its AIC to that of all the others
aics = AIC(l.exp,l.weib,l.lognorm,l.gamma,l.loglogis)
aics[order(aics$AIC),]




##############################################################################
##############################################################################

# HW 12

p = read.csv("/Users/jamescutler/Desktop/Survival/pbcirrhosis.csv")

# DON'T FORGET TO CHANGE PBC'S STATUS VARIABLE LEVEL 1 to 0 and 2 TO 1
p$status[p$status == 1] = 0
p$status[p$status == 2] = 1

p$edemaFact = as.character(p$edema)
summary(p$edemaFact)
p$edemaFact = as.factor(p$edemaFact)
summary(p$edemaFact)

p.weib = survreg(Surv(futime,status) ~ edemaFact, data = p, dist = "weibull")
summary(p.weib)
-2*-1711
p.exp = survreg(Surv(futime,status) ~ edemaFact, data = p, dist = "exponential")
summary(p.exp)

cont.weib = survreg(Surv(futime,status) ~ edema, data = p, dist = "weibull")
cont.exp = survreg(Surv(futime,status) ~ edema, data = p, dist = "exponential")
anova(cont.weib,cont.exp)

myLRT(p.weib$loglik[2],p.exp$loglik[2],1,2)
anova(p.exp,p.weib)

# Weibull is the better model (we can reject the null of accepting the nested).


# # relevel just for fun (to see if it's the same as SAS)
# p.rev =
#   p %>%
#   mutate(edemaFact = relevel(edemaFact, ref = "1"))
# p.rev$edemaFact
# rev.weib = survreg(Surv(futime,status) ~ edemaFact, data = p.rev, dist = "weibull")
# rev.exp = survreg(Surv(futime,status) ~ edemaFact, data = p.rev, dist = "exponential")
# summary(rev.weib)



# Interpret the effect of edema from the Weibull model
summary(p.weib)
e0.5medR = exp(-0.8118)
e1.0medR = exp(-2.0506)
1-e0.5medR
1-e1.0medR

# IF YOU EXPONENTIATE THE ESTIMATE ("Value") IN THE MODEL, YOU GET THE MEDIAN SURVIVAL
## TIME RATIO. IF YOU EXPONENTIATE THE NEGATIVE OF THE ESTIMATE, YOU GET THE HAZARD RATIO.
e0.5HR = exp(0.8118/.864)
e1.0HR = exp(2.0506/.864)
e0.5HR
e1.0HR


# Compare edema continuous Weibull with edema categorical weibull 
## YOU CAN'T COMPARE THEM USING LRT! - THEY'RE NOT NESTED! Source: https://stats.stackexchange.com/questions/181337/why-is-this-nested
## You can compare them using AIC, BIC, etc.
## WAIT A MINUTE!!!! ACTUALLY, DR. DING SAYS THE CONTINUOUS **IS** NESTED WITHIN THE 
## CATEGORICAL!!!!!! His explanation makes sense too. The categorical will have more
## coefficients, and it could be linear, or it could be bent in some way (non-linear).
## So if we treat it as a straight line in a linear model (i.e., make it continuous),
## then we are REDUCING it.
myLRT(cont.weib$loglik[2],p.weib$loglik[2],2,2)
anova(cont.weib,p.weib)

# It looks like we can accept the continuous model (P=0.28).
summary(cont.weib)
eMedR = exp(-1.9575)
eHRscale = exp(1.9575/0.868)
eHR = exp(1.9575)
1-eMedR
eHRscale
eHR

pdt = data.table(p)
pdt[edema == 0,median(futime)]
pdt[edema == .5,median(futime)]
pdt[edema == 1,median(futime)]
1-299/1882
1 - 1180.5/1882



# Cox Snell residual plot
## Source: http://biostat.mc.vanderbilt.edu/wiki/pub/Main/QingxiaChen/Lab12.txt

## edema coefficient
mu.p = cont.weib$coef[1]
## scale
sigma.p = cont.weib$scale
## Hazard ratio for edema 1 compared to edema 0
lambda.p = exp(-mu.p/sigma.p)
## 1/scale
alpha.p = 1/sigma.p
## log of hazard ratio
beta.p = -cont.weib$coef[2]/sigma.p
## CS residuals and plot
cs.p = lambda.p*exp(with(p, cbind(edema))%*%beta.p)*(p$futime^alpha.p)
s.cs.p = survfit(Surv(cs.p, p$status) ~ 1, type = "fleming-harrington")
plot(s.cs.p, mark.time = FALSE, fun = "cumhaz",
     xlab = "Cox-Snell Residuals", ylab = "Estimated Cumulative Hazard Rates",
     main = "Cox-Snell Residual Plot for Continuous Edema",
     conf.int = FALSE)
abline(a=0,b=1, col = "blue")
