

# Note: A lot of this code is not clean either. It's sometimes pretty hard to tell
## why I'm doing this or that. Sorry!


# Analysis of Frequency Data



# GREAT SOURCES OF INFORMATION
## The Stats Geek
### On Hosmer-Lemeshow GOF for logistic regression: https://thestatsgeek.com/2014/02/16/the-hosmer-lemeshow-goodness-of-fit-test-for-logistic-regression/
## R-Bloggers
### On evaluating logistic regression models: https://www.r-bloggers.com/evaluating-logistic-regression-models/
## R Companion for the Handbook of Biological Statistics, Salvatore Mangiafico
### On logistic regression - EXCELLENT SOURCE: https://rcompanion.org/rcompanion/e_06.html





# Libraries
library(readxl)
library(sas7bdat)
library(ggplot2)
library(tidyr)
library(broom)
library(data.table)

# library(ggstatsplot)   # Mother of all packages - plots only work in base R
library(lubridate)
library(scales)
library(wrapr)           # For build_frame

# library(emmeans)       # Which one is better, emmeans or lme4?
# library(lsmeans)
# library(lme4)

# THE MUST-HAVE EPI/FREQUENCY PACKAGES
# Believe it or not, the mantelhaen.test (Cochran-Mantel-Haenszel Chi-Squared Test For Count Data) is in stats - base R!
library(DescTools)       # For CochranArmitageTest; BreslowDayTest
library(epitools)        # duh
library(epiR)            # For sample size calculation for CC studies
library(rawr)            # Use rawr's jt.test instead of clinfun's
library(lazyWeave)       # For Mantel-Haenszel test
library(questionr)       # For chi square standardized residuals - VERY IMPORTANT!
# THE MUST-HAVE LOGISTIC REGRESSION PACKAGES
library(car)                    # For VIF (variance inflation factor) to detect collinearity
## and for Wald test chi-square statistic?
library(ResourceSelection)      # For Hosmer-Lemeshow: hoslem.test()
library(pROC)                   # For ROC (AUC, or C statistic)
library(LogisticDx)             # For overall C stat and ROC, Pearson, and deviance

# library(clinfun)       # install.packages("clinfun"); for jonckheere.test - NO GOOD - use rawr's jt.test instead
library(CATT)            # install.packages("CATT"); cochran-armitage test







# Custom functions
## Goodness of linear fit test
## DAT file string to data frame (first n columns can have character values)
## Almighty RR and OR (witch CIs and P-values) marginal table generator (makes these ^ obsolete)
## AR, AR%, PAR, PAR%
## Format regression model output into Estimate, CIs, P-value dataframe
## pvalPretty
## Plot stacked barchart proportions for CC study strata and crude
## countsToDF
## Overlapping barplot for visualizing trends from ix2 tables
## Mantel-Haenszel test with row scores and column scores (cuz lazyWeave's didn't work)





# Chi-square goodness of linear fit function
## Requires: library(lazyWeave)
goflinear = function(contingencyTable){
  ct.chis = chisq.test(contingencyTable)
  ct.mant = mantel.test(contingencyTable)
  new.chis = ct.chis$statistic - ct.mant$statistic
  gof.pvalue = 1 - pchisq(new.chis,df=ct.chis$parameter-1)
  print("And the chi-square goodness of linear fit p-value is ...")
  return(gof.pvalue)
}


# String data to dataframe, with first n columns having character values
myDATstring.to.df = function(mystring,mystringNames,n=0){
  mystring = unlist(strsplit(mystring,"\n"))
  mystring = gsub("^\\s+ ","",mystring)
  mystring = gsub("\\s+$","",mystring)
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
  for (i in (1+n):ncol(myDF)){
    myDF[,i] = as.numeric(as.character(myDF[,i]))
  } 
  return(myDF)
}


# Convert a Desc(tabs) object to a full-on marginal table with all the extras
myDescToAlmightyDF = function(Desc2b2b2,OutcomeCol1,OutcomeCol2,Confounder){
  print("Must have yesOC first. Must have yesR first.")
  myPercent = 3
  m2b2b2 = capture.output(Desc2b2b2)
  nums = m2b2b2[12:17]
  nums = gsub("\\s+"," ",nums)
  mylength = sapply(strsplit(nums," "),length)      # genius code
  numlines = vector()
  for (i in 1:length(nums)){
    numlines[i] = 
      paste(strsplit(nums," ")[[i]][(mylength[i]-2):(mylength[i])],collapse = " ")
  }
  # numlines = gsub("\\s+"," ",numlines)
  numlines = unlist(strsplit(numlines," "))
  a = numlines[seq(1,length(numlines),3)]     # Will become outcome column 1
  b = numlines[seq(2,length(numlines),3)]     # Will become outcome column 2
  c = numlines[seq(3,length(numlines),3)]     # Will become sum column
  numdf = data.frame(OutcomeCol1 = as.numeric(as.character(a)),
                     OutcomeCol2 = as.numeric(as.character(b)),
                     "Sum" = as.numeric(as.character(c)))
  # Compute row percents
  numdf$rowPercent = round((numdf$OutcomeCol1 / numdf$Sum)*100, myPercent)
  # Compute column percents for outcome status 1
  numdf$colPercent1 = NA
  numdf$colPercent1[1] = round( numdf[1,1]/sum(numdf[1,1],numdf[2,1]) ,myPercent)
  numdf$colPercent1[3] = round( numdf[3,1]/sum(numdf[3,1],numdf[4,1]) ,myPercent)
  numdf$colPercent1[5] = round( numdf[5,1]/sum(numdf[5,1],numdf[6,1]) ,myPercent)
  # Compute column percents for outcome status 2
  numdf$colPercent2 = NA
  numdf$colPercent2[1] = round( numdf[1,2]/sum(numdf[1,2],numdf[2,2]) ,myPercent)
  numdf$colPercent2[3] = round( numdf[3,2]/sum(numdf[3,2],numdf[4,2]) ,myPercent)
  numdf$colPercent2[5] = round( numdf[5,2]/sum(numdf[5,2],numdf[6,2]) ,myPercent)
  # Exposure column
  expos = c(strsplit(nums," ")[[1]][2],strsplit(nums," ")[[2]][2])    # This line is unreliable--does the second strsplit need to be a [[2]][2] or a [[2]][1]?
  # Confounder column
  confound = c(strsplit(nums," ")[[1]][1],strsplit(nums," ")[[3]][1])
  # Data frame of the marginal table
  martabDF = expand.grid(Exposure = expos,
                         Confounder = confound)
  martabDF = martabDF[,2:1]
  martabDF[,1] = as.character(martabDF[,1])
  ## add sum rows to df
  martabDF=
    rbind(martabDF,
          c("Sum",expos[1]),
          c("Sum",expos[2]))
  ## add outcome columns to df
  martabDF$OutcomeCol1 = numdf$OutcomeCol1
  martabDF$OutcomeCol2 = numdf$OutcomeCol2
  ## add sum column to df
  martabDF$Sum = numdf$Sum
  ## add row percents column (for Cohort) to df
  martabDF$rowPercent = numdf$rowPercent
  ## add column percents columns (for CC) to df
  martabDF$colPercent1 = numdf$colPercent1
  martabDF$colPercent2 = numdf$colPercent2
  ## make the column names good
  colnames(martabDF) = c(Confounder,"Exposure",OutcomeCol1,OutcomeCol2,"Sum",
                         "Cohort: rowPct","CC: P(E|D)","CC: P(E|ND)")
  
  # From marginal table to 3 individ tables with measures of assoc., their CIs, and P-values
  allrownames = martabDF$Exposure[1:2]
  allcolnames = colnames(martabDF)[3:4]
  # Partial table 1
  tab1 = vector()
  for (i in 1:2){
    tab1[i] = martabDF[1,i+2]
  }
  for (i in 3:4){
    tab1[i] = martabDF[2,i]
  }
  tab1 = matrix(tab1,2,2,byrow = TRUE)
  dimnames(tab1) = list("Exposure" = allrownames, "Outcome" = allcolnames)
  # Partial table 2
  tab2 = vector()
  for (i in 1:2){
    tab2[i] = martabDF[3,i+2]
  }
  for (i in 3:4){
    tab2[i] = martabDF[4,i]
  }
  tab2 = matrix(tab2,2,2,byrow = TRUE)
  dimnames(tab2) = list("Exposure" = allrownames, "Outcome" = allcolnames)
  # Crude (sum) table
  tab3 = vector()
  for (i in 1:2){
    tab3[i] = martabDF[5,i+2]
  }
  for (i in 3:4){
    tab3[i] = martabDF[6,i]
  }
  tab3 = matrix(tab3,2,2,byrow = TRUE)
  dimnames(tab3) = list("Exposure" = allrownames, "Outcome" = allcolnames)
  
  # So far we have martabDF, and tab1, tab2, and tab3
  # Now to get RRs and ORs for each tab (1,2,3), their CIs and P-values, and add them to martabDF
  
  # ORs
  or.tab1 = oddsratio.wald(tab1,rev = "both")
  or.tab2 = oddsratio.wald(tab2,rev = "both")
  or.tab3 = oddsratio.wald(tab3,rev = "both")
  # RRs
  rr.tab1 = riskratio.wald(tab1,rev = "both")
  rr.tab2 = riskratio.wald(tab2,rev = "both")
  rr.tab3 = riskratio.wald(tab3,rev = "both")
  
  # Add an Exposure OR column for CC, with CI columns and p.value
  martabDF$CC_E.OR = NA
  martabDF$CC_E.OR[1] = or.tab1$measure[2,1]      # The odds ratio for partial table 1
  martabDF$CC_E.OR[3] = or.tab2$measure[2,1]      # The odds ratio for partial table 2
  martabDF$CC_E.OR[5] = or.tab3$measure[2,1]      # The odds ratio for the crude table
  martabDF$E.OR_CI.low = NA
  martabDF$E.OR_CI.low[1] = or.tab1$measure[2,2]
  martabDF$E.OR_CI.low[3] = or.tab2$measure[2,2]
  martabDF$E.OR_CI.low[5] = or.tab3$measure[2,2]
  martabDF$E.OR_CI.upp = NA
  martabDF$E.OR_CI.upp[1] = or.tab1$measure[2,3]
  martabDF$E.OR_CI.upp[3] = or.tab2$measure[2,3]
  martabDF$E.OR_CI.upp[5] = or.tab3$measure[2,3]
  martabDF$E.OR_p.value = NA
  martabDF$E.OR_p.value[1] = or.tab1$p.value[2,3]
  martabDF$E.OR_p.value[3] = or.tab2$p.value[2,3]
  martabDF$E.OR_p.value[5] = or.tab3$p.value[2,3]
  # Add a Disease RR column for Cohort, with CI columns and p.value
  martabDF$RR = NA
  martabDF$RR[1] = rr.tab1$measure[2,1]           # The risk ratio for partial table 1 
  martabDF$RR[3] = rr.tab2$measure[2,1]           # The risk ratio for partial table 2 
  martabDF$RR[5] = rr.tab3$measure[2,1]           # The risk ratio for the crude table 
  martabDF$RR_CI.low = NA
  martabDF$RR_CI.low[1] = rr.tab1$measure[2,2]
  martabDF$RR_CI.low[3] = rr.tab2$measure[2,2]
  martabDF$RR_CI.low[5] = rr.tab3$measure[2,2]
  martabDF$RR_CI.upp = NA
  martabDF$RR_CI.upp[1] = rr.tab1$measure[2,3]
  martabDF$RR_CI.upp[3] = rr.tab2$measure[2,3]
  martabDF$RR_CI.upp[5] = rr.tab3$measure[2,3]
  martabDF$RR_p.value = NA
  martabDF$RR_p.value[1] = rr.tab1$p.value[2,3]
  martabDF$RR_p.value[3] = rr.tab2$p.value[2,3]
  martabDF$RR_p.value[5] = rr.tab3$p.value[2,3]
  
  # Now round everything
  z = c(1,3,5)
  assoc.round = 4
  for (i in z){
    martabDF$CC_E.OR[i] = round(martabDF$CC_E.OR[i],assoc.round)
    martabDF$E.OR_CI.low[i] = round(martabDF$E.OR_CI.low[i],assoc.round)
    martabDF$E.OR_CI.upp[i] = round(martabDF$E.OR_CI.upp[i],assoc.round)
    martabDF$E.OR_p.value[i] = round(martabDF$E.OR_p.value[i],assoc.round)
    
    martabDF$RR[i] = round(martabDF$RR[i],assoc.round)
    martabDF$RR_CI.low[i] = round(martabDF$RR_CI.low[i],assoc.round)
    martabDF$RR_CI.upp[i] = round(martabDF$RR_CI.upp[i],assoc.round)
    martabDF$RR_p.value[i] = round(martabDF$RR_p.value[i],assoc.round)
  }
  ## -c(6,13:16) if you want a CC study table
  ## -c(7:12)    if you want a Cohort study table
  print("Select -c(6,13:16) if you want a CC study table")
  print("Select -c(7:12)    if you want a Cohort study table")
  return(martabDF)
}


# Risk difference (attributable risk) and AR%, PAR, and PAR%
myARPAR = function(my2x2,abcd.vec=NULL){
  print("Requires epitools")
  if (is.null(abcd.vec)){
    my.prop = prop.table(my2x2,1)
    my.mar = table.margins(my2x2)
    # AR = Ie - Iue
    Ie = my.prop[1,1]
    Iue = my.prop[2,1]
    ar = Ie - Iue
    # AR% = AR/Ie
    arpct = ar/Ie
    # PAR = Ipop - Iue 
    Ipop = my.mar[3,1]/my.mar[3,3]
    par = Ipop - Iue
    # PAR% = PAR/Ipop
    parpct = par/Ipop
  } else{
    d = c("D","ND")
    e = c("E","UE")
    de = matrix(abcd.vec,2,2,byrow = TRUE)
    dimnames(de) = list("Exposure" = e, "Disease" = d)
    de.pr = prop.table(de,1)
    de.mar = table.margins(de)
    # AR = Ie - Iue
    Ie = de.pr[1,1]
    Iue = de.pr[2,1]
    ar = Ie - Iue
    # AR% = AR/Ie
    arpct = ar/Ie
    # PAR = Ipop - Iue
    Ipop = de.mar[3,1]/de.mar[3,3]
    par = Ipop - Iue
    # PAR% = PAR/Ipop
    parpct = par/Ipop
  }
  # WHY DOESN'T THIS WORK?
  # z = c(ar,arpct,par,parpct)
  # for (i in z){
  #   z[i] = round(z[i],3)
  # }
  ar = round(ar,3)
  arpct = round(arpct,3)
  par = round(par,3)
  parpct = round(parpct,3)
  # Interpretations: AR, AR%, PAR, PAR%
  print(sprintf("%s cases of [disease] per 1,000 persons with the exposure can be attributed to their exposure.",ar*1000))
  print(sprintf("%s percent of [disease] among those with the exposure is due to their exposure.",arpct*100))
  print(sprintf("%s cases of [disease] per 1,000 persons in the population can be attributed to the exposure.",par*1000))
  print(sprintf("%s percent of [disease] in the population is due to the exposure.",parpct*100))
  print("AR, AR%, PAR, and PAR%:")
  return(c(ar,arpct,par,parpct))
}


# Format the results of your regression model in a nice convenient way (Estimate, CIs, P-value)
logis.mod.results = function(yourmodel,CI=.95){
  print("Confidence intervals are Wald confidence intervals (same as SAS).")
  ORs = exp(yourmodel$coefficients[-1])
  allCIs = exp(confint.default(yourmodel,level = CI))
  CIs.low = allCIs[2:(length(allCIs)/2)]
  CIs.high = allCIs[(length(allCIs)/2 + 2):length(allCIs)]
  tmod = tidy(yourmodel)
  Pvals = tmod$p.value[-1]
  myresults = data.frame(OR = ORs,
                         CI.low = CIs.low,
                         CI.high = CIs.high,
                         P.value = Pvals)
  return(myresults)
}


# Convert p-values to pretty p-values
pvalPretty = function(pvals){
  prettyVals = prettyNum(pvals,scientific = FALSE)
  return(matrix(prettyVals))
}


# Plot proportions P(E|D) and P(E|ND) each stratum at a time, and crude as well
## ONLY WORKS FOR CC STUDIES - EXPOSURE IS THE FILL AND X GROUPING IS DISEASE STATUS
CCProportions = function(my2x2,
                         values = c("#F8766D","#00BFC4"),
                         lab.Ex.Pos.first,
                         perc.round = 3,
                         mytitle,
                         mysubtitle,
                         subtitle.color){
  my2x2trans = t(my2x2)
  my2x2melt = melt(my2x2trans,
                   varnames = c("Disease","Exposure"),
                   id.vars = "Disease")
  ggplot(my2x2melt %>% group_by(Disease) %>%                       # CC groups by disease ...
           mutate(Percent = round(value/sum(value),perc.round)),
         aes(x = Disease, y = Percent,                             # ... thus CC has disease on x ...
             fill = Exposure,cumulative = TRUE)) +                 # ... and CC has exposure as fill
    geom_col() +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = values,
                      labels = lab.Ex.Pos.first) +
    geom_text(aes(label = paste0(Percent*100,"%")),
              position = position_stack(vjust = .5)) +
    theme(plot.subtitle = element_text(color = subtitle.color)) +
    labs(title = mytitle,
         subtitle = mysubtitle,
         x = "", y = "Percentage")
}


# For going from counts data frames to de-tabled dataframes (long dataframes)
countsToCases <- function(df, countcol = "Freq") {
  # Get the row indices to pull from df
  idx <- rep.int(seq_len(nrow(df)), df[[countcol]])
  
  # Drop count column
  df[[countcol]] <- NULL
  
  # Get the rows from df
  df[idx, ]
}


# For creating a REAL FREAKING data frame out of a table
## Designed for a kappa2() type situation. You plug the dataframe that results 
## from this function into the kappa2() function.
tbToDF = function(atable,dim1names,dim2names){
  tbtovec = as.vector(t(atable))
  dim2rep = rep(dim2names,nrow(atable))
  mtable = table.margins(atable)
  rater1.mars = mtable[1:(nrow(mtable)-1),ncol(mtable)]
  rater1 = capture.output(
    for(i in 1:nrow(atable)){
      cat(rep(dim1names[i],rater1.mars[i])," ")
    }
  )
  rater1 = unlist(strsplit(rater1," "))
  rater1 = rater1[-which(rater1 == "")]
  
  # Now for rater 2
  rater2 = capture.output(
    for (i in 1:length(tbtovec)){
      cat(rep(dim2rep[i],tbtovec[i])," ")
    }
  )
  rater2 = unlist(strsplit(rater2," "))
  rater2 = rater2[-which(rater2 == "")]
  df = as.data.frame(cbind(rater1,rater2))
  return(df)
}

# Overlapping barplot for visualizing trends from ix2 tables
overlap.bars = function(df,dfcol1,dfcol2,xbreakLabs=1:nrow(df),title,xlab,ylab,
                        colors = c("green","blue"),subtitle=NULL){
  ggplot() +
    geom_bar(data = df, mapping = aes(1:nrow(df),dfcol1),
             stat = "identity", width = .9, fill = colors[1]) +
    geom_bar(data = df, mapping = aes(1:nrow(df),dfcol2),
             stat = "identity", width = .4, fill = colors[2]) +
    scale_x_continuous(breaks = 1:nrow(df),
                       labels = xbreakLabs) +
    labs(title = title, subtitle = subtitle,
         x=xlab,y=ylab)
}


# Mantel-Haenszel test with row.scores and col.scores (because lazyWeave's didn't work)
myMantel.haenszel = function(temp,rowscores,colscores){
  u <- rowscores
  v <- colscores
  n <- sum(temp)
  p1 <- sum(t(t(temp) * v) * u)
  p2 <- sum(u * rowSums(temp))
  p3 <- sum(v * colSums(temp))
  p4 <- sum(u^2 * rowSums(temp))
  p5 <- sum(u * rowSums(temp))^2
  p6 <- sum(v^2 * colSums(temp))
  p7 <- sum(v * colSums(temp))^2
  num <- p1 - (p2 * p3/n)
  den <- sqrt((p4 - p5/n) * (p6 - p7/n))
  r <- num/den
  M2 <- (n - 1) * r^2
  pvalue <- 1 - stats::pchisq(M2, 1)
  print("M^2 and p-value")
  return(c(M2,pvalue))
}

# Getting (rough) midpoints for interval row levels
midpoints = function(df.interval.row.levels){
  print("Note! This requires intervals with a dash in between the two termini.")
  mychar = as.character(df.interval.row.levels)
  mychar = unlist(strsplit(mychar,"-"))
  mynum = as.numeric(mychar)
  mids = vector()
  l = length(mynum)/2
  for (i in 1:l){
    mids[i] = median(c(mynum[2*i],mynum[(2*i)-1]))
  }
  return(mids)
}



###############################################################################
###############################################################################
###############################################################################

### END CUSTOM FUNCTIONS ###

###############################################################################
###############################################################################
###############################################################################
















###############################################################################
###############################################################################
###############################################################################

### SART ###

###############################################################################
###############################################################################
###############################################################################

# Trying to figure out teacher salary regression

data("mtcars")
?mtcars
colnames(mtcars)
summary(mtcars)

sal = read.csv("/Users/jamescutler/Desktop/Data_Course_cutler/Salaries_Professors.csv")
sal = sal[,-1]
summary(sal)
sal$discipline = as.character(sal$discipline)
sal$discipline[which(sal$discipline == "A")] = "theory"
sal$discipline[which(sal$discipline == "B")] = "applied"
sal$discipline = as.factor(sal$discipline)
summary(sal)


mod.gen = glm(salary ~ rank+discipline+yrs.since.phd+yrs.service+sex,
              data = sal)
summary(mod.gen)

sex.means = lsmeans(mod.gen, "sex")
plot(sex.means)
plot(sal$yrs.service ~ sal$sex)
plot(sal$yrs.since.phd ~ sal$sex)

service.lm = lm(salary ~ yrs.service, data = sal)
plot(sal$yrs.service,sal$salary, col = sal$sex)
abline(service.lm)

sincePhD.lm = lm(salary ~ yrs.since.phd, data = sal)
plot(sal$yrs.since.phd,sal$salary, col = sal$sex)
abline(sincePhD.lm)

summary(sincePhD.lm)
summary(service.lm)

mod.lm = lm(salary ~ yrs.service+yrs.since.phd+sex, 
            data = sal)
summary(mod.lm)
sex.yrs.means = lsmeans(mod.lm, "sex")
plot(sex.yrs.means)
sex.sal.means = tapply(sal$salary,sal$bin.sex,mean)

sal$bin.sex = NA
sal$bin.sex[which(sal$sex == "Male")] = 1
sal$bin.sex[which(sal$sex == "Female")] = 0

basic.sex.aov.mod = aov(salary ~ sex, 
                        data = sal)
summary(basic.sex.aov.mod)

ggplot(sal, aes(x=sex,y=salary)) +
  geom_violin()


plot(sal$bin.sex,sal$salary)
points(x=c(0,1),y=c(sex.sal.means[1][[1]],sex.sal.means[2][[1]]), 
       pch = 18, col = c("green","red"))
abline(h=c(sex.sal.means[1][[1]],sex.sal.means[2][[1]]),
       col = c("green","red"))
abline(h=c(106080.2,114537.2),col = c("darkgreen","purple"))
sex.sal.means
sex.means
qqnorm(sal$salary)
qqline(sal$salary, col = "red")

gen.means = plot(sex.means) + 
  coord_cartesian(xlim = c(88000,118000))
yrs.means = plot(sex.yrs.means)
grid.arrange(gen.means,yrs.means)





#############################################################################
#############################################################################

# Day 1 (8/19) exercise in class

(90/70)/(403/1201)
(90*1201)/(70*403)

# The hard way to do it:
hc = data.frame(htn = rep(c("present","absent"),c(493,1271)),
                cd = 
                  c(
                    rep(c("present","absent"),c(90,403)),
                    rep(c("present","absent"),c(70,1201)) 
                  ) 
)
hc.t = table(hc$htn,hc$cd)

# The easy way to do it:
htn = c("highBP","normalBP")
cvd = c("hasCVD","noCVD")
htn.cvd = matrix(c(90,403,
                   70,1201),2,2,byrow = TRUE)
dimnames(htn.cvd) = list("HTN exposure" = htn, "CVD" = cvd)
htn.cvd

# Measures of association:
rr = riskratio(htn.cvd, rev = "both")
oddsratio(htn.cvd, rev = "both")
oddsratio.fisher(htn.cvd, rev = "both")
or = oddsratio.wald(htn.cvd, rev = "both")    # THIS IS THE RIGHT ONE!

# The proportions of the useful one will show you CI_exp and CI_unexp:
prop.table(htn.cvd,2)
prop.table(htn.cvd,1) # This is the useful one
.18255/.05507         # This is the RR

rr$measure            # Same as this
or$measure            # Why should the OR be higher than the RR?
or$measure[2,3]

# Now to look at the conf. ints. graphically:
rr.vs.or = function(or.measure,rr.measure,exposure,disease){
  or = or.measure[2,1]
  or.lb = or.measure[2,2]
  or.ub = or.measure[2,3]
  rr = rr.measure[2,1]
  rr.lb = rr.measure[2,2]
  rr.ub = rr.measure[2,3]
  df = data.frame(Ratio = c("RR","OR"),
                  vals = c(rr,or),
                  lb = c(rr.lb,or.lb),
                  ub = c(rr.ub,or.ub))
  ggplot(df, aes(x=Ratio,y=vals,col=Ratio)) +
    geom_point(aes(fill=Ratio),size=4,shape=21) +
    geom_errorbar(aes(ymin=lb,ymax=ub),width=.3) +
    theme_bw() +
    labs(x="Type of ratio",y="",
         title = sprintf("OR is %s and RR is %s",round(or,4),round(rr,4)),
         subtitle = sprintf("Association between %s and %s",exposure,disease))
}
rr.vs.or(or$measure,rr$measure,"hypertension","CVD") # Perfect.




##########################################################

# HW 1:

# The funny thing about problem 1.3 is that I didn't even remember that I should
## know the mean and variance of any binomial random variable as long as I know
## the number of trials and the probability of success, because it's np and 
## np(1-p) duh!
n = 100
p = .25
vr = n*p*(1-p)
vr
n*p
s.d = sqrt(vr)
s.d
d = dbinom(x=50,size=100,prob=.25)
dbinom(x=25,size=100,prob=.25)
test = rbinom(10000,size=100,prob=.25)
hist(test)
pbinom(q=50,size=100,prob=.25,lower.tail = FALSE)
qbinom(p=d,size=100,prob=.25,lower.tail = FALSE)
# This is weird stuff. Is the probability 2.13e-8 or 4.5e-8???

# Problem 1.8:
yes = 344
n = 1170
yes/n
prop.test(x=yes,n=n,p=.5, alternative = "less")

BinomCI(x=yes,n=n, conf.level = .99, method = "clopper-pearson")
CI.for1prop = function(CI,prop1,n1){
  alpha = 1 - CI
  z = qnorm(1 - alpha/2)
  me = z*sqrt( prop1*(1-prop1) / n1 )
  ub = prop1 + me
  lb = prop1 - me
  cat(sprintf("We are %s percent confident that the true population proportion is \nbetween %.4f and %.4f",CI*100,lb,ub))
}
CI.for1prop(.99,yes/n,n)

m = c("wilson", "wald", "agresti-coull", "jeffreys",
      "modified wilson", "wilsoncc","modified jeffreys",
      "clopper-pearson", "arcsine", "logit", "witting", "pratt")
for (i in 1:length(m)) {
  print(m[i])
  print(BinomCI(x=yes,n=n, conf.level = .99, method = m[i]))
}
# It looks like the Wald is the one I can replicate.

# last problem:
htn = c("highBP","normalBP")
cvd = c("hasCVD","noCVD")
htn.cvd = matrix(c(90,403,
                   70,1201),2,2,byrow = TRUE)
dimnames(htn.cvd) = list("HTN exposure" = htn, "CVD" = cvd)
htn.cvd


hi = c("injured","notInjured")
wh = c("yes","noHelmet")
wh.hi = matrix(c(17,218,
                 130,428),2,2,byrow=TRUE)
dimnames(wh.hi) = list("Head injury" = hi, "Wore helmet" = wh)
wh.hi
235/793
646/793
(235 + 646 - 218)/793
218/646
chisq.test(wh.hi, correct = FALSE)

pi = 235/793
q = 1-pi
pi/q
iGy = 17/147
iGy/(1-iGy)
iGn = 218/646
iGn/(1-iGn)



#########################
rf = c("E","UE")
out = c("D","ND")
rfout = matrix(c(60,40,
                 68,32),2,2,byrow=TRUE)
dimnames(rfout) = list("Smoking" = rf, "Cancer" = out)
rfout

chisq.test(rfout)
chisq.test(rfout,correct = FALSE)
or = oddsratio.wald(rfout, rev = "both")
or$measure
rr = riskratio.wald(rfout, rev="both")
rr$measure
riskratio(rfout, rev = "both")
riskratio.small(rfout, rev="both")
riskratio.boot(rfout, rev="both")
odd

prop.table(rfout,2)
prop.table(rfout,1)

rfout





#############################################################################

# The probability of getting any number of heads out of 
## 2,5,10,20,100,1000 coin tosses:

two.coin = dbinom(x=0:2,size=2,prob=.5)
plot(two.coin, type = "l",lwd=2,col="red", ylim=c(0,1))
abline(h=c(0,.05), col = "brown")

five.coin = dbinom(x=0:5,size=5,prob=.5)
plot(five.coin, type = "l",lwd=2,col="orange", ylim=c(0,1))
abline(h=c(0,.05), col = "brown")

ten.coin = dbinom(x=0:10,size=10,prob=.5)
plot(ten.coin, type = "l",lwd=2,col="yellow", ylim=c(0,1))
abline(h=c(0,.05), col = "brown")

twenty.coin = dbinom(x=0:20,size=20,prob=.5)
plot(twenty.coin, type = "l",lwd=2,col="green", ylim=c(0,1))
abline(h=c(0,.05), col = "brown")

hundo.coin = dbinom(x=0:100,size=100,prob=.5)
plot(hundo.coin, type = "l",lwd=2,col="blue", ylim=c(0,1))
abline(h=c(0,.05), col = "brown")
hundo.coin[80]
sum(hundo.coin[1:100])
or.better = sum(hundo.coin[70:100])
prettyNum(1/or.better,big.mark = ",")
# See Notes app STATS QUESTIONS to see my questions about this.

thous.coin = dbinom(x=0:1000,size = 1000,prob = .5)
h.thou = head(tail(thous.coin,510),20)
plot(h.thou, col = "purple", ylim = c(0,.05))
plot(thous.coin, type = "l",lwd=2,col="purple", ylim=c(0,1))
abline(h=c(0,.05), col = "brown")
# This one teaches you that the chance of getting a certain number of heads isn't all there 
## is to deciding whether the event is surprising or not. You should also check how the 
## probability of the event compares to any other possible events. See below:
length(thous.coin) # IT'S 1001! NOT 1000!
mymax = which(thous.coin == max(thous.coin))
thous.coin[mymax:1001] # It's symmetric, so just take one half (indices 501 to 1001).
thou.df = data.frame(probs = thous.coin[1001:mymax],
                     percentiles = seq(0,100,length.out = length(thous.coin[1001:mymax]))
)
plot(thou.df$percentiles,thou.df$probs, cex = .5)
abline(v=c(90,92,94),col=c("red","orange","green"))
comfortable.prob.range = thou.df$probs[which(thou.df$percentiles == 94)]
which(thous.coin >= comfortable.prob.range)
less.comfort.range = thou.df$probs[which(thou.df$percentiles == 92)]
which(thous.coin >= less.comfort.range)
1/sum(thous.coin[1:460])
# That's what I'm talking about.

# Show them all together now (it doesn't line up because of the varying index):
plot(thous.coin, type = "l",lwd=2,col="purple", ylim=c(0,1))
points(two.coin, type = "l",lwd=2,col="red", ylim=c(0,1))
points(five.coin, type = "l",lwd=2,col="orange", ylim=c(0,1))
points(ten.coin, type = "l",lwd=2,col="yellow", ylim=c(0,1))
points(twenty.coin, type = "l",lwd=2,col="green", ylim=c(0,1))
points(hundo.coin, type = "l",lwd=2,col="blue", ylim=c(0,1))





#############################################################################
#############################################################################

62.4-1.3
62.4/1.3
prettyNum(61.1/1e6,scientific = FALSE)

0.001304+0.000121

o1 = 0.001304/(1-0.001304)
o2 = 0.000121/(1-0.000121)
o1/o2



inc = c("low.inc","high.inc")
coverage = c("yes coverage","no coverage")
inc.cov = matrix(c(1617,678,
                   3317,356),2,2,byrow=TRUE)
dimnames(inc.cov) = list("Income" = inc, "Coverage" = coverage)
inc.cov


prop.table(inc.cov,1)
prop.table(inc.cov,2)

or = oddsratio.wald(inc.cov, rev = "both")
or$measure


CI.for1prop = function(CI,prop1,n1){
  alpha = 1 - CI
  z = qnorm(1 - alpha/2)
  me = z*sqrt( prop1*(1-prop1) / n1 )
  ub = prop1 + me
  lb = prop1 - me
  print(sprintf("We are %s percent confident that the true population proportion is between %.4f and %.4f",CI*100,lb,ub))
}
CI.for1prop(.95,0.7045752,3673)  # WHY ARE THEY SLIGHTLY DIFFERENT?  
prop.test(1617,1617+678)         # WHY ARE THEY SLIGHTLY DIFFERENT?

or$p.value
rr = riskratio(inc.cov, rev = "both")
rr$measure

or$p.value
prop.table(inc.cov,1)
prop.table(inc.cov,2)

prop.test(inc.cov)
0.7045752 - 0.9030765






#############################################################################
#############################################################################

# Week 3 in class

party = c("Democrat","Republican","Independent")
gender = c("Women","Men")
pg = matrix(c(495,272,590,
              330,265,498),nrow=2,ncol=3,byrow=TRUE)
dimnames(pg) = list("Gender"=gender,"Party" = party) # Define row dim FIRST, col second
pg
chisq.test(pg)
chisq.test(pg[,c(1,3)])
prop.table(pg[,c(1,3)],1)
prop.table(pg,1)
DI = rowSums(pg[,c(1,3)])
R = pg[,2]
DI
R

chisq.test(x=DI,y=R)
class(DI)
DI.R = cbind(DI,R)
chisq.test(DI.R)
prop.table(DI.R,1)


# Sample size with epiR's epi.ccsize function:








###############

# Week 3 guided exercise

cancer = c("Controlled","Not Contr.")
tx = c("Surgery","Radiation")
cancer.tx = matrix(c(21,2,
                     15,3),2,2,byrow = TRUE)
dimnames(cancer.tx) = list("Treatment" = tx,"Cancer" = cancer)
cancer.tx
fisher.test(cancer.tx)

party = c("Democrat","Republican","Independent")
race = c("White","Black")
pr = matrix(c(871,821,336,
              347,42,83),2,3,byrow = TRUE)
dimnames(pr) = list("Race" = race, "Party" = party)
pr
prop.table(pr,1)
prop.table(pr,2)



# Team exercise

# 1a
epi.ccsize(OR = 3, 
           p0 = .15, 
           n = NA, 
           power = .8, 
           r = 1,
           conf.level = .95)
# 1b
epi.ccsize(OR = 3, 
           p0 = .10, 
           n = NA, 
           power = .8, 
           r = 1,
           conf.level = .95)

# 2a
epi.ccsize(OR = 2.5, 
           p0 = .25, 
           n = NA, 
           power = .8, 
           r = 1,
           conf.level = .95)

# 2b
epi.ccsize(OR = 2.5, 
           p0 = .35, 
           n = NA, 
           power = .8, 
           r = 1,
           conf.level = .95)

# 3a
epi.ccsize(OR = 2,
           p0 = .25,
           n = NA,
           power = .9)
# 3b
epi.ccsize(OR = 3,
           p0 = .25,
           n = NA,
           power = .9)
# 3c
epi.ccsize(OR = 2,
           p0 = .16,
           n = NA,
           power = .9)
# 3d
epi.ccsize(OR = 3,
           p0 = .16,
           n = NA,
           power = .9)

# 4a
epi.propsize(treat=.6,
             control = .5,
             n = NA,
             power = .95,
             conf.level = .99)

# 5a
epi.propsize(treat = .2,
             control = .05,
             n = NA,
             power = .8,
             conf.level = .9)
# 5b
epi.propsize(treat = .2,
             control = .05,
             n = NA,
             power = .8,
             conf.level = .95)

# 6a - COULDN'T FIGURE THIS ONE OUT
epi.ccsize(OR = NA,
           p0 = .1,
           n = 100000002,
           power = .8,
           r = 2/99999999)








#############################################################################
#############################################################################



# HW 3 a
a = 
  epi.propsize(treat = .2,
               control = .35,
               n = NA,
               power = .9,
               conf.level = .95)
b = 
  epi.propsize(treat = .2,
               control = .4,
               n = NA,
               power = .9,
               conf.level = .95)
c = 
  epi.propsize(treat = .2,
               control = .45,
               n = NA,
               power = .9,
               conf.level = .95)
a$n.total
b$n.total
c$n.total




######################

# Week 3 Report

# SEE Titanic_logistic_ORs




####################################################################################
####################################################################################

# Mantel-Haenszel chi-square test for linear trend
# Cochran-Armitage trend test
# Jonckheere-Terpstra test
# A chi-square goodness of linear fit test

# Mantel-Haenszel first:
alcohol = c("0","<1","1-2","3-5","≥6")
malform = c("No","Yes")
alcmal = matrix(c(17066,48,
                  14464,38,
                  788,5,
                  126,1,
                  37,1),5,2,byrow = TRUE)
dimnames(alcmal) = list("Alcohol" = alcohol, "Malformation" = malform)
alcmal
mantel.test(alcmal)
prop.table(alcmal,1)
.0263/.0028

CochranArmitageTest(alcmal, alternative = "increasing")
?mantel.test
?mantelhaen.test
?CochranArmitageTest
?prop.trend.test
?prop.test
prop.test(alcmal, alternative = "greater")
mantelhaen.test(alcmal)

# To do prop.trend.test with alcmal:
as.matrix(rowSums(alcmal))
alcmal
n = rowSums(alcmal)
x = alcmal[,2]
x
n
prop.trend.test(x,n)


# What do they think Jonckheere-Terpstra is?
set.seed(1234)
g = rep(1:5, rep(10,5))
x = rnorm(50)
g
x
plot(1:50,x+.8*g)
?JonckheereTerpstraTest()

wtf = Untable(c(4,6,5), type = "ordered")[,1]
wtf
coffee = data.frame(time = c(447,396,383,410,438,521,468,391,
                             504,472,513,543,506,489,407),
                    group = wtf)
coffee
plot(coffee$time ~ coffee$group)
JonckheereTerpstraTest(time ~ group, data = coffee)



# Team exercise:

tb = "0 (<25) 	59 	133 
1 (26-35) 	54 	150 
2 (36-45)	53 	120 
3 (46-55) 	61 	120 
4 (>55) 	61 	57 
"
tb = unlist(strsplit(tb,"\n"))
tb = unlist(strsplit(tb,"\t"))
ages = tb[seq(1,length(tb),3)]
ages = gsub(" ","",ages)

Yes = tb[seq(2,length(tb),3)]
No = tb[seq(3,length(tb),3)]
No = gsub(" ","",No)
Yes = gsub(" ","",Yes)
No = as.numeric(No)
Yes = as.numeric(Yes)
No
Yes


TB = c("Yes","No")
age = ages
tbage = cbind(Yes,No)
dimnames(tbage) = list("Age" = age, "TB" = TB)
tbage
prop.table(tbage,1)
prop.table(tbage,1)[,1]
barplot(prop.table(tbage,1)[,1])

chisq.test(tbage)
mantel.test(tbage)
CochranArmitageTest(tbage, alternative = "decreasing")

# For prop.trend.test
n = rowSums(tbage)
x = tbage[,1]
x
prop.trend.test(x,n) # This is useless. It doesn't mean it's nonlinear. The alternative is that there is a linear trend.

# Now try with a legitimate test for linearity:
tbl=matrix(c(20,10,20,20,20,30),nrow=2)
tbl
CATT(table = tbl)
# CATT requires a 2xJ table, not an Ix2 table:
CATT(table = t(tbage)) # I don't think this works. Why are the two cochran tests giving me complete different values?
t(tbage)




#####################################################

# HW 4

inc = "(<5, 5–15, 15–25,>25)"
inc = unlist(strsplit(inc,","))
inc = gsub("\\(","",inc)
inc = gsub("\\)","",inc)
inc = gsub(" ","",inc)
inc
tally = "(2, 4, 13, 3 / 2, 6, 22, 4 / 0, 1, 15, 8 / 0, 3, 13, 8)"
satis = c("very dis","little satis","mod satis","very satis")
tally = unlist(strsplit(tally,","))
tally = unlist(strsplit(tally,"/"))
tally = gsub("\\(","",tally)
tally = gsub("\\)","",tally)
tally = gsub(" ","",tally)
tally = as.numeric(tally)
job = matrix(tally,nrow=4,ncol=4,byrow = TRUE)
dimnames(job) = list("Income" = inc, "Satisfaction" = satis)
job

chisq.test(job)
std.resids = tidy(chisq.residuals(job,std = TRUE))        # chisq.residuals is from questionr
# setwd("/Users/jamescutler/Desktop/Frequency_Data/")
# write.csv(std.resids,"standardized_residuals.csv")

rowscores = c(3,10,20,35)
colscores = c(1,3,4,5)
mantel.test(job,
            row.scores = rowscores,
            col.scores = colscores)
lazyWeave:::mantel.test
nrow(job)
ncol(job)

data("mtcars")
mtcars$gear
mtcars$cyl
table(mtcars$gear,mtcars$cyl)
mantel.test(x=mtcars$gear,byVar = mtcars$cyl,
            row.scores = c(1,2,3),
            col.scores = c(1,2,3))
temp = job


myMantel.haenszel(job,rowscores = rowscores,colscores = colscores)


# Vet study
vet = "45-49      1      1      0      4
45-49      1      1      1      1
45-49      1      0      0      17
45-49      1      0      1      6
50-54      2      1      0      3
50-54      2      1      1      7
50-54      2      0      0      29
50-54      2      0      1      11
55-59      3      1      0      3
55-59      3      1      1      9
55-59      3      0      0      15
55-59      3      0      1       6
60-64      4      1      0      2
60-64      4      1      1      10
60-64      4      0      0      12
60-64      4      0      1      5
65-69      5      1      0      10
65-69      5      1      1      4
65-69      5      0      0      6
65-69      5      0      1      10
70-74      6      1      0      3
70-74      6      1      1      3
70-74      6      0      0      3
70-74      6      0      1      3"
# cols = c("Age","Age Group","Case","Exposed","Counts")
cols = "Age Age_Group Case Exposed Counts"

vet = myDATstring.to.df(vet,cols,n=1)
vetdt = data.table(vet)
smoke = vetdt[,sum(Counts),by=c("Age","Exposed")]
smoke

smoke.tab = 
  spread(data = smoke,
         key = Exposed,
         value = V1)
smoke.tab
colnames(smoke.tab) = c("Age","NonSmoker","Smoker")
smoke.tab

# Now to visualize the trends:
overlap.bars(smoke.tab,smoke.tab$NonSmoker,smoke.tab$Smoker,xbreakLabs = smoke.tab$Age,
             title = "Smokers (blue) and Non-Smokers (green) by Age",
             xlab = "Age Group", ylab = "Counts in each exposure status category")

# For case status:
case = vetdt[,sum(Counts),by=c("Age","Case")]
case.tab = 
  spread(case,
         key = Case,
         value = V1)
case.tab
overlap.bars(case.tab,case.tab$`0`,case.tab$`1`,
             case.tab$Age,title = "Case status (blue=case,green=control) by age",
             xlab="Age Group",ylab = "Counts in each case status category")

# Tests:
CochranArmitageTest()        # a frequency table or a matrix, alternative = t,i,d
jonckheere.test()            # x and g, alternative = t,i,d
JonckheereTerpstraTest()     # x and g, alternative = t,i,d
mantel.test()                # x, equal/midrank
myMantel.haenszel()          # table, rowscores,colscores

# Now convert smoke.tab to a contingency table:
s.tab = 
  cbind(smoke.tab$`0`,smoke.tab$`1`)
dimnames(s.tab) = list("Age" = smoke.tab$Age,
                       "Smoking Status" = c("No","Yes"))
s.tab
CochranArmitageTest(s.tab, alternative = "two.sided")   # Gives double the (right?) one
CochranArmitageTest(s.tab, alternative = "increasing")  # Gives the right one?
CochranArmitageTest(s.tab, alternative = "d")           # Gives the compliment

mantel.test(s.tab)                                      # Gives nearly the same as two.sided Cochrane
myrowscores = midpoints(smoke.tab$Age)
myMantel.haenszel(s.tab,
                  rowscores = c(45,50,55,60,65,70),
                  colscores = c(1,2))                   # Same as mantel.test!
myMantel.haenszel(s.tab,
                  rowscores = myrowscores,
                  colscores = c(1,2))                   # Same as mantel.test! 


# The strange SAS-universe concept of a "chi-square goodness of linear fit" test:
## Coffee dataset example from class:
class.coffee = "0	yes	20
1-2	yes	153
3-4	yes	106
>=5	yes	88
0	no	88
1-2	no	271
3-4	no	154
>=5	no	130
"
classC = 
  charCol1.text.to.df(class.coffee,"Cups Disease Counts",n=2)
c.tab = 
  spread(classC,
         key = Disease,
         value = Counts)
c.tab = c.tab[c(2,3,4,1),]
c.tab
overlap.bars(c.tab,c.tab$no,c.tab$yes,xbreakLabs = c.tab$Cups,
             title = "Disease status (green=no disease, blue=disease) by coffe consumption",
             x="Daily cups of coffee",y="Disease status counts")

# c.tab (Coffee dataset) as a contingency table:
cctab = 
  cbind(c.tab$no,c.tab$yes)
dimnames(cctab) = list("Cups" = c.tab$Cups, "Disease" = c("Yes","No"))
cctab

# Now to get the goodness of linear fit p-value:
cctab.chis = chisq.test(cctab)
cctab.mant = mantel.test(cctab)
new.chis = cctab.chis$statistic - cctab.mant$statistic
# AND NOW BEHOLD I HAVE MY GOODNESS OF LINEAR FIT P-VALUE:
gof.pvalue = 
  1-pchisq(new.chis,df=cctab.chis$parameter-1)
# Now to create a function for it:
gof(cctab)


# Now use the gof function on the vet study data:
smoke.tab
s.tab
gof(s.tab)


# Agresti 2.19:
rel = "Less that High School 178 (4.5) 138 (−2.6) 108 (−1.9)
High School or Junior College 570 (2.6) 648 (1.3) 442 (−4.0)
Bachelor or Graduate 138 (−6.8) 252 (0.7) 252 (6.3)"
rel = unlist(strsplit(rel,"\n"))
rel.cats = gsub("([0-9]+).*$","",rel)
rel.cats = gsub(" $","",rel.cats)
rel.cats
# gsub("[^0-9]","",rel)   # Weird.
rel.nums = str_remove(rel,rel.cats)
rel.nums = gsub("\\([^\\)]+\\)","",rel.nums)
rel.nums = gsub("\\s+"," ",rel.nums)
rel.nums = gsub("^ ","",rel.nums)
rel.nums = gsub(" $","",rel.nums)
rel.nums = as.numeric(unlist(strsplit(rel.nums," ")))
rel.nums

reltab = matrix(rel.nums,3,3,byrow = TRUE)
dimnames(reltab) = list("Education" = rel.cats, 
                        "Religiosity" = c("Fundamentalist","Moderate","Liberal"))
reltab

tidy(chisq.residuals(reltab,std = TRUE))
chisq.test(reltab)

# JT test:
reltab
epitools::table.margins(reltab)                # TABLE MARGINS!
jonckheere.test(mtcars$gear,mtcars$cyl)

# install.packages("remotes")
# remotes::install_github("raredd/rawr")
# library(rawr)
jt.test(reltab)  # HECK YEAH!!!! Same z statistic as SAS










####################################################################################
####################################################################################

# Week 5 - three-way contingency tables (IxJxK tables)
## Using x=array(data.vector,IxJxK), dimnames(x), names(dimnames(x)) to make 3-way tables
## Using xtabs and expand.grid to make 3-way tables (BETTER THAN ARRAY)
## Cochran-Mantel-Haenszel chi-squared Test for conditional independence 
### Source: https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/mantelhaen.test
## Breslow-Day test for homogeneity of odds ratios
### Source: https://rdrr.io/cran/DescTools/man/BreslowDayTest.html
## Using xtabs(), cbind() and expand.grid() with Desc() and my DescToDF() for 3-way table!

pol = read.csv("/Users/jamescutler/Desktop/Frequency_Data/PolIdeolData.csv")
pol

?table
?array
array(pol,c(2,2,5))

array(1:20,c(2,2,5))
tab = array(1:8,c(2,2,2))
dimnames(tab) = list(c('No','Yes'),c('No','Yes'),c('No','Yes'))
names(dimnames(tab)) = c("Exposure","Disease","Confounder")
tab

smoking = c(42,165,7,28,97,65,436,474)

s = array(smoking,c(2,2,2))
dimnames(s) = list(c('Yes','No'),c('Yes','No'),c('Yes','No'))
names(dimnames(s)) = c("Exposure","Disease","Age ≥ 65")
s
smar = table.margins(s)
smar
smar[1]



# Simpson's paradox (confounding) in comparison of kidney stone treatments
tx = c("open","PNL")
outcome = c("success","failure")
# Group 1
one = matrix(c(81,6,
               234,36),2,2,byrow = TRUE)
dimnames(one) = list("Treatment" = tx,
                     "Outcome" = outcome)
one
or1 = oddsratio.wald(one,rev="both")
or1$p.value
or1$measure

# Group 2
two = matrix(c(192,71,
               55,25),2,2,byrow=TRUE)
dimnames(two) = list("Treatment" = tx,
                     "Outcome" = outcome)
two
or2 = oddsratio.wald(two,rev="both")
or2$p.value
or2$measure

# Totals
tots = matrix(c(273,77,
                289,61),2,2,byrow = TRUE)
dimnames(tots) = list("Treatment" = tx,
                      "Outcome" = outcome)
tots
ortot = oddsratio.wald(tots,rev="both")      # Holy cow!! Simpson's paradox!!! wow!!!!!!
ortot$p.value
ortot$measure




# DescTools' BreslowDayTest for homogeneity of odds ratios
## Source: https://rdrr.io/cran/DescTools/man/BreslowDayTest.html - excellent example copied below!
# Questions:
## 1) Does it have to be a 2x2xK table for Breslow-Day to work? If so, why?
migraine <- xtabs(freq ~ .,
                  cbind(expand.grid(treatment=c("active", "placebo"),
                                    response =c("better", "same"),
                                    gender   =c("female", "male")),
                        freq=c(16, 5, 11, 20, 12, 7, 16, 19))
)
?xtabs

# Amazing:
expand.grid(treatment=c("active", "placebo"),
            response =c("better", "same"),
            gender   =c("female", "male"))
# Amazing:
cbind(expand.grid(treatment=c("active", "placebo"),
                  response =c("better", "same"),
                  gender   =c("female", "male")),
      freq=c(16, 5, 11, 20, 12, 7, 16, 19))

# EVEN MORE AMAZING! A 3-WAY TABLE!!
migraine

# And the result:
BreslowDayTest(migraine)


# Now to do the smoking 3-way table with xtabs and expand.grid!
## The trick is to get the order of the columns right, based on what you see 
## in the lecture (slide 21). This order is different, so the freq vector
## order has to be different.
expand.grid(Exposure = c("Smoker","Non-Smoker"),
            Outcome = c("Dead","Not Dead"),
            Age = c("≥ 65","< 65"))
cbind(
  expand.grid(Exposure = c("Smoker","Non-Smoker"),
              Outcome = c("Dead","Not Dead"),
              Age = c("≥ 65","< 65")),
  freq = c(42,165,7,28,97,65,436,474)
)

smoke = 
  xtabs(freq ~ .,
        cbind(
          expand.grid(Exposure = c("Smoker","Non-Smoker"),
                      Outcome = c("Dead","Not Dead"),
                      Age = c("≥ 65","< 65")
          ),
          freq = c(42,165,7,28,97,65,436,474)
        )
  )
smoke
Desc(smoke)
# This ^ tells me that for Desc, it goes first row var, second row var, column var
# So let's try it again with the variables in the same order as the lecture:
expand.grid(Age = c("< 65","≥ 65"),
            Exposure = c("Smoker","Non-Smoker"),
            Outcome = c("Dead","Not Dead"))
freq = c(97,42,65,165,436,7,474,28)
cbind(
  expand.grid(Age = c("< 65","≥ 65"),
              Exposure = c("Smoker","Non-Smoker"),
              Outcome = c("Dead","Not Dead")),
  freq = c(97,42,65,165,436,7,474,28)
)

# xtabs is the best! xtabs, expand.grid() and cbind()
smoke = 
  xtabs(freq ~ .,
        cbind(
          expand.grid(Age = c("<65","≥65"),
                      Exposure = c("Smoker","Non-Smoker"),
                      Outcome = c("Dead","Not Dead")
          ),
          freq = c(97,42,65,165,436,7,474,28)
        )
  )
smoke
Desc(smoke)

smoke.conf.tab = 
  myDescToDF(Desc(smoke),"Dead","Not.Dead","Age") 
smoke.conf.tab                # Hell. Yeah.


smoke.individ = confDFtoTables(smoke.conf.tab)
smoke.individ[[1]]
for (i in 1:3){
  print(riskratio(smoke.individ[[i]],rev="both")$measure) 
}
# rev="both" is absolutely necessary to get the right answer here!!

# So, the conditional RRs are 1.51 and 1.00, and the marginal RR is .76, so ...
## ... naturally, the weighted overall RR should be something in that ballpark:
mantelhaen.test(smoke)

# Try rearranging the tables
smoke.partial.tabs = margin.table(smoke,c(2,3,1))
mantelhaen.test(smoke.partial.tabs) # BTW, the continuity correction only changes the p-value, not the OR
smoke.partabs2 = margin.table(smoke,c(3,2,1))
mantelhaen.test(smoke.partabs2)

# Try the cmh.test from lawstat
install.packages("lawstat")
library(lawstat)
smoke.cmh=
  cmh.test(smoke.partial.tabs)
smoke.cmh$parameter

# Try EpiR
epi.2by2(dat = smoke, method = "cohort.count")
epi.2by2(dat = smoke.partial.tabs, method = "cohort.count")
riskratio.wald(smoke)

# Try out the titanic dataset
data("Titanic", package = "datasets")
str(Titanic)
Titanic
# indices: 1=class; 2=gender; 3=age; 4=mortality
partial.tabs = margin.table(Titanic, c(2,4,1))   # So it goes c(row,column,table#)
# The c(2,4,1) ensures that it's a 2x2xK table, where K is the only dimension that 
## can be greater than 2. So class has to be the last dimension, otherwise you would
## have a 4x2x2 or a 2x4x2. In this case, it has to be 2x2x4, so only the sex and
## survival variables can alternate positions.
partial.tabs
mantelhaen.test(partial.tabs)
apply(partial.tabs,3,oddsratio.wald)









###############################################################################
###############################################################################

# Report 2 - CHD in males with catecholamine levels and age covariates

chd = read_excel("/Users/jamescutler/Desktop/Frequency_Data/CHDCAT.xlsx")
chd$CHD[chd$CHD == "PRESENT"] = "CHD"
chd$CHD[chd$CHD == "ABSENT"] = "No CHD"
chd$CAT[chd$CAT == "PRESENT"] = "High_CAT"  # Has to be one word
chd$CAT[chd$CAT == "ABSENT"] = "Low_CAT"    # Has to be one word
chd

# Make INDIVIDUAL partial tables (SEE plots below)

# Make BOTH PARTIAL TABLES TOGETHER
chd.tabs = 
  xtabs(N ~ .,
        chd)
chd.tabs

# Now make the marginal table out of the both together partial tables
Desc(chd.tabs)


# MAKE A PUBLICATION-READY DATAFRAME OUT OF IT
newCHDMarTab = 
  myDescToAlmightyDF(Desc(chd.tabs),"CHD","No_CHD","Age")
newCHDMarTab

# Select the right columns for a CC study
## -c(6,13:16) if you want a CC study table
## -c(7:12)    if you want a Cohort study table
newCHD.martab.final = 
  newCHDMarTab[,-c(6,13:16)]

# Polish the column names if necessary
colnames(newCHD.martab.final) = c("Age","Exposure","CHD","No CHD","Sum",
                                  "% High CAT | CHD","% High CAT | No CHD",
                                  "OR","95% CI lower","95% CI upper","P-value")

# Polish confounder and exposure level labels if necessary
newCHD.martab.final$Exposure = gsub("\\_"," ",newCHD.martab.final$Exposure)
newCHD.martab.final$Age = c("≥ 55","≥ 55",
                            "< 55","< 55",
                            "Sum","Sum")

# Ready for kable
newCHD.martab.final
save(newCHD.martab.final, file = "/Users/jamescutler/Desktop/Frequency_Data/CHD_marginal_table")




# The concept behind the Breslow-Day test for interaction makes sense. An interaction
## is where the effect of one variable (here, CAT) changes depending upon the status
## of the other variable (here, age group). 
## Ho: All stratum-specific ORs are equal.
## Ha: Not all stratum-specific ORs are equal.

# In order to get the right Breslow-Day statistic andCMH common (weighted) OR for 
## this CC, you have to make sure you split the tables by age, rather than by outcome
chd.tabs
chd.cmh.tabs = margin.table(chd.tabs,c(2,3,1)) # Now age, which was the row var, is the table index
chd.cmh.tabs
BreslowDayTest(chd.cmh.tabs)    # Same as SAS Breslow-Day p-value (P=0.0069).
# There is an interaction; not all stratum-specific ORs are equal.

# CMH analysis
mantelhaen.test(chd.cmh.tabs) # 1.804 matches the adjusted CMH OR from SAS (CI: 1.17-2.78)




# Plots!
barplot(prop.table(chd.age.mar,2))
barplot(prop.table(chd.lt55.tab,2))
barplot(prop.table(chd.ge55.tab,2))

# The real barplots
library(reshape2)
# From https://stackoverflow.com/questions/46894576/ggplot-grouped-barplot-for-2x2-table
data <- matrix(c(4, 1, 5, 2), ncol = 2, dimnames = list(c("C", "D"), c("A", "B")))
data_m <- melt(data, varnames = c("Exp", "Obs"), id.vars = "Exp")

ggplot(data_m %>% group_by(Exp) %>% 
         mutate(perc = round(value/sum(value),2)), 
       aes(x = Exp, y = perc, 
           fill = Obs, cumulative = TRUE)) +
  geom_col() +
  geom_text(aes(label = paste0(perc*100,"%")), 
            position = position_stack(vjust = 0.5))




# Partial tables for the plots of each stratum
chd.lt55 = chd[chd$AGE=="LT55",]
chd.ge55 = chd[chd$AGE=="GE55",]
# Partial table 1
chd.lt55.tab = xtabs(N ~ CAT+CHD,
                     data=chd.lt55)
chd.lt55.tab
# Partial table 2
chd.ge55.tab = xtabs(N ~ CAT+CHD,
                     data = chd.ge55)
chd.ge55.tab

# Adjustable variables to plug into all graphs
mytitle = "Percent high CAT levels by disease status"
values = c("#F8766D","#00BFC4")
labels = c("High CAT","Low CAT")
perc.round = 3

# Less than 55 stratum
lt55.plot = t(chd.lt55.tab)

chdlt55Melt = melt(lt55.plot,
                   varnames = c("Disease","Exposure"),
                   id.vars = "Disease")

# 55 and greater stratum
ge55.plot = t(chd.ge55.tab)

ge55Melt = melt(ge55.plot,
                varnames = c("Disease","Exposure"),
                id.vars = "Disease")


# Plot
stratum.lt55=
  ggplot(chdlt55Melt %>% group_by(Disease) %>% 
           mutate(Percent = round(value/sum(value),perc.round)),
         aes(x = Disease, y = Percent,
             fill = Exposure,cumulative = TRUE)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = values,
                    labels = labels) +
  geom_text(aes(label = paste0(Percent*100,"%")),
            position = position_stack(vjust = .5)) +
  theme(plot.subtitle = element_text(color = "purple")) +
  labs(title = mytitle,
       subtitle = "For those less than 55-yo",
       x = "", y = "Percentage")
stratum.lt55

# Plot
stratum.ge55=
  ggplot(ge55Melt %>% group_by(Disease) %>% 
           mutate(Percent = round(value/sum(value),perc.round)),
         aes(x = Disease, y = Percent,
             fill = Exposure,cumulative = TRUE)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = values,
                    labels = labels) +
  geom_text(aes(label = paste0(Percent*100,"%")),
            position = position_stack(vjust = .5)) +
  theme(plot.subtitle = element_text(color = "purple")) +
  labs(title = mytitle,
       subtitle = "For those 55-yo and greater",
       x = "", y = "Percentage")


# Marginal
chd.age = margin.table(chd.tabs,c(2,3,1))
chd.age.mar = margin.table(chd.age,c(1,2))
chdSum.plot = t(chd.age.mar)

chdSumMelt = melt(chdSum.plot,
                  varnames = c("Disease","Exposure"),
                  id.vars = "Disease")
# Plot
chd.marginal=
  ggplot(chdSumMelt %>% group_by(Disease) %>% 
           mutate(Percent = round(value/sum(value),perc.round)),
         aes(x = Disease, y = Percent,
             fill = Exposure,cumulative = TRUE)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = values,
                    labels = labels) +
  geom_text(aes(label = paste0(Percent*100,"%")),
            position = position_stack(vjust = .5)) +
  theme(plot.subtitle = element_text(color = "red")) +
  labs(title = mytitle,
       subtitle = "Marginal percentages",
       x = "", y = "Percentage")
chd.marginal

# Beautiful
grid.arrange(stratum.lt55,
             stratum.ge55,
             chd.marginal,
             ncol=3)







###############################################################################
###############################################################################

# Risk difference (attributable risk) and AR%, PAR, PAR%

# See myARPAR function at top under list of custom functions

d = c("CVD","No CVD")
e = c("HBP","No HBP")
de = matrix(c(90,403,
              70,1201),2,2,byrow = TRUE)
dimnames(de) = list("Exposure" = e,"Disease" = d)
de

# This function is so awesome
myARPAR(de)
myARPAR(abcd.vec = c(90,403,
                     70,1201))









###############################################################################
###############################################################################

# Review for Midterm Exam

make2x2fromVector = function(abcd.vec){
  d = c("D","ND")
  e = c("E","UE")
  de = matrix(abcd.vec,2,2,byrow = TRUE)
  dimnames(de) = list("Exposure" = e, "Disease" = d)
  
  return(de) # Can't return this in table form if it's inside a combine-into-a-vector function
}
de = 
  make2x2fromVector(c(71,72,
                      72,75))
table.margins(de)

# Doing it by hand
e11 = 143*143/290
e12 = 147*143/290
e21 = 143*147/290
e22 = 147*147/290
(71-e11)^2 /e11 +
  (72-e12)^2 /e12 +
  (72-e21)^2 /e21 +
  (75-e22)^2 /e22
1-pchisq(.01304757,1)







###############################################################################
###############################################################################

# Week 8 Logistic Regression!

# Seatbelt dataset
car = 
  "1 1 1   34
1 1 0   273
1 0 1   87
1 0 0   221
0 1 1    6
0 1 0   37
0 0 1   12
0 0 0   51"

car = 
  myDATstring.to.df(car,"seat belt injury n")
car
car$seat[car$seat == 1] = "Front"
car$seat[car$seat == 0] = "Back"
car$belt[car$belt == 1] = "With"
car$belt[car$belt == 0] = "Without"
car$injury[car$injury == 1] = "Injured"
car$injury[car$injury == 0] = "Safe"
car

car.tabs = xtabs(n ~ .,
                 data = car)
car.tabs

car.desc = 
  Desc(car.tabs)

car.mar = 
  myDescToAlmightyDF(car.desc,"Injured","Safe","Seat")

# CC car.mar
car.cc = car.mar[,-c(6,13:16)]

# Cohort car.mar
car.cohort = car.mar[,-c(7:12)]
car.cohort  

# CMH adjusted OR
car.cmh.tabs = margin.table(car.tabs,c(2,3,1))
car.cmh.tabs
mantelhaen.test(car.cmh.tabs)


# Logistic regression
car.mod1 

# CAD dataset
cad = 
  "1	1	1	21
1	1	0	6
1	0	1	9
1	0	0	9
0	1	1	8
0	1	0	10
0	0	1	4
0	0	0	11"
cad = myDATstring.to.df(cad,"sex EKG CAD n")
cad$sex[cad$sex == 1] = "Male"
cad$sex[cad$sex == 0] = "Female"
cad$EKG[cad$EKG == 1] = "Depression"
cad$EKG[cad$EKG == 0] = "No_Depress"
cad$CAD[cad$CAD == 1] = "CAD"
cad$CAD[cad$CAD == 0] = "No_CAD"
cad

cad.tabs = xtabs(n ~ .,
                 data = cad)
cad.tabs

cad.desc = Desc(cad.tabs)
cad.desc

cad.mar = myDescToAlmightyDF(cad.desc,"CAD","No_CAD","Sex")
cad.mar

# CC cad.mar
cad.cc = cad.mar[,-c(6,13:16)]
cad.cc

# Cohort cad.mar
cad.cohort = cad.mar[,-c(7:12)]
cad.cohort

# CAD - CMH adjusted OR
cad.cmh.tabs = margin.table(cad.tabs,c(2,3,1))
cad.cmh.tabs
mantelhaen.test(cad.cmh.tabs)


# Logistic regression
library(MASS)
data("menarche")
menarche
cbind(menarche$Menarche,menarche$Total - menarche$Menarche)





# In class exercise
hos = 
  "Hospital A 	New Treatment 	Improved 		93 
Hospital A 	New Treatment 	Did Not Improve 	147 

Hospital A 	Standard Treatment 	Improved 		59 
Hospital A 	Standard Treatment 	Did Not Improve 	181 

Hospital B 	New Treatment 	Improved 		41 
Hospital B 	New Treatment 	Did Not Improve 	79 

Hospital B 	Standard Treatment 	Improved 		33 
Hospital B 	Standard Treatment 	Did Not Improve 	87 
"
hos = gsub("Hospital ","",hos)
hos = gsub(" Treatment","",hos)
cat(hos)

hos = 
  "A 1 1 93
A	1 0 147
A 0 1 59
A 0 0 181
B	1 1 41
B	1 0 79
B	0 1 33
B	0 0 87"
hos = gsub("\t"," ",hos)
hos 
cat(hos)
hos = myDATstring.to.df(hos,"hosptial tx effect n",n=1)
hos


hos$tx[hos$tx == 1] = "New"
hos$tx[hos$tx == 0] = "Standard"
hos$effect[hos$effect == 1] = "Improved"
hos$effect[hos$effect == 0] = "Not_Improved"
hos

hos.tabs = xtabs(n ~.,
                 hos)
hos.tabs

hos.desc = Desc(hos.tabs)
hos.desc

hos.mar = myDescToAlmightyDF(hos.desc,"Improved","Not_Improved","Hospital")
hos.cc = hos.mar[,-c(6,13:16)]
hos.cohort = hos.mar[,-c(7:12)]
hos.cc

# Interaction and confounding
hos.cmh.tabs = margin.table(hos.tabs,c(2,3,1))
hos.cmh.tabs
BreslowDayTest(hos.cmh.tabs)
mantelhaen.test(hos.cmh.tabs)


# Logistic regression
hos.df = as.data.frame(hos.tabs)

hos.mo1 = glm(effect ~ tx + hosptial, weights = Freq, data = hos.df,
              family = binomial(link = "logit"))
summary(hos.mo1)

# ORs for hospital and tx
exp(.03941) # For hospital
exp(.54654) # For tx

coefs = hos.mo1$coefficients[-1]
class(coefs)
dimnames(coefs)
hos.mo1$terms[[1]][3]
exp(confint(hos.mo1))

log.mod.results = function(yourmodel,CI=.95){
  ORs = exp(yourmodel$coefficients[-1])
  CIs = exp(confint(yourmodel,level = CI))[-1,]
  Pvals 
}


exp(confint(hos.mo1,level = .95))[-1,]


tidy(hos.mo1)

or.ci = 
  ggcoefstats(
    x = hos.mo1,
    exponentiate = TRUE,
    output = "tidy"
  )
or.ci[order(or.ci$estimate, decreasing = TRUE),]
or.ci


or.ci[,c(1,2,3,4,9)]








###################################################################################
###################################################################################

# Week 8 HW

d2 = 
  "1 1 1   98
1 1 0   102
1 0 1   53
1 0 0   147
2 1 1   11
2 1 0   49
2 0 1   12
2 0 0   48"
d2 = myDATstring.to.df(d2,"Hospital  TX  Improve n")
d2

d3 =
  "1 1 1   81
1 1 0   54
1 0 1   31
1 0 0    9
2 1 1    9
2 1 0   36
2 0 1   40
2 0 0   93"
d3 = myDATstring.to.df(d3,"Hospital  TX  Improve n")
d3


hos = 
  "A 1 1 93
A	1 0 147
A 0 1 59
A 0 0 181
B	1 1 41
B	1 0 79
B	0 1 33
B	0 0 87"
hos = gsub("\t"," ",hos)
hos 
cat(hos)
hos = myDATstring.to.df(hos,"hosptial tx effect n",n=1)
hos

d1 = hos
colnames(d1) = colnames(d2)
d1$Hospital = as.character(d1$Hospital)
d1$Hospital[d1$Hospital == "A"] = 1
d1$Hospital[d1$Hospital == "B"] = 2
d3
d2
d1


# Improve rename
d1$Improve[d1$Improve == 1] = "Improved"
d1$Improve[d1$Improve == 0] = "Not_Improved"

d2$Improve[d2$Improve == 1] = "Improved"
d2$Improve[d2$Improve == 0] = "Not_Improved"

d3$Improve[d3$Improve == 1] = "Improved"
d3$Improve[d3$Improve == 0] = "Not_Improved"

# Treatment rename
d1$TX[d1$TX == 1] = "New"
d1$TX[d1$TX == 0] = "Standard"

d2$TX[d2$TX == 1] = "New"
d2$TX[d2$TX == 0] = "Standard"

d3$TX[d3$TX == 1] = "New"
d3$TX[d3$TX == 0] = "Standard"




# tabs and Descs and marginal tables

d1.tabs = xtabs(n ~ .,
                data = d1)
d1.tabs

d2.tabs = xtabs(n ~ .,
                data = d2)
d2.tabs

d3.tabs = xtabs(n ~ .,
                data = d3)
d3.tabs

d1.desc = Desc(d1.tabs)
d2.desc = Desc(d2.tabs)
d3.desc = Desc(d3.tabs)

d1.mar = myDescToAlmightyDF(d1.desc,"Improved","Not_Improved","Hospital")
d1.cc = d1.mar[,-c(6,13:16)]

d2.mar = myDescToAlmightyDF(d2.desc,"Improved","Not_Improved","Hospital")
d2.cc = d2.mar[,-c(6,13:16)]

d3.mar = myDescToAlmightyDF(d3.desc,"Improved","Not_Improved","Hospital")
d3.cc = d3.mar[,-c(6,13:16)]

d1.cc
d2.cc
d3.cc




# Breslow and CMH 
d1.cmh.tabs = margin.table(d1.tabs,c(2,3,1))
d2.cmh.tabs = margin.table(d2.tabs,c(2,3,1))
d3.cmh.tabs = margin.table(d3.tabs,c(2,3,1))

BreslowDayTest(d1.cmh.tabs)
BreslowDayTest(d2.cmh.tabs)
BreslowDayTest(d3.cmh.tabs)

mantelhaen.test(d1.cmh.tabs)
mantelhaen.test(d2.cmh.tabs)
mantelhaen.test(d3.cmh.tabs)




# Logistic
d1.df = as.data.frame(d1.tabs)
d1.df
d1.mod = glm(Improve ~ TX + Hospital, weights = Freq, data = d1.df,
             family = binomial(link = "logit"))
summary(d1.mod)

d2.df = as.data.frame(d2.tabs)
d2.df
d2.mod = glm(Improve ~ TX + Hospital, weights = Freq, data = d2.df,
             family = binomial(link = "logit"))
summary(d2.mod)

d3.df = as.data.frame(d3.tabs)
d3.df
d3.mod = glm(Improve ~ TX + Hospital, weights = Freq, data = d3.df,
             family = binomial(link = "logit"))
summary(d3.mod)

tidy(d3.mod)

# The ultimate model results formatter
logis.mod.results = function(yourmodel,CI=.95){
  ORs = exp(yourmodel$coefficients[-1])
  allCIs = exp(confint(yourmodel,level = CI))
  CIs.low = allCIs[2:(length(allCIs)/2)]
  CIs.high = allCIs[(length(allCIs)/2 + 2):length(allCIs)]
  tmod = tidy(yourmodel)
  Pvals = tmod$p.value[-1]
  myresults = data.frame(OR = ORs,
                         CI.low = CIs.low,
                         CI.high = CIs.high,
                         P.value = Pvals)
  return(myresults)
}


d1.modResults = logis.mod.results(d1.mod)

d2.modResults = logis.mod.results(d2.mod)

d3.modResults = logis.mod.results(d3.mod)

d1.modResults
d2.modResults
d3.modResults




# Relevel the variables so the interpretation is correct
d1.rev.df = 
  d1.df %>%
  mutate(Improve = relevel(Improve, ref = "Not_Improved")) %>%
  mutate(TX      = relevel(TX,      ref = "Standard")) %>%
  mutate(Hospital= relevel(Hospital,ref = 2))

d2.rev.df = 
  d2.df %>%
  mutate(Improve = relevel(Improve, ref = "Not_Improved")) %>%
  mutate(TX      = relevel(TX,      ref = "Standard")) %>%
  mutate(Hospital= relevel(Hospital,ref = 2))

d3.rev.df = 
  d3.df %>%
  mutate(Improve = relevel(Improve, ref = "Not_Improved")) %>%
  mutate(TX      = relevel(TX,      ref = "Standard")) %>%
  mutate(Hospital= relevel(Hospital,ref = 2))




# Correctly-reference-grouped models
d1.modrev = glm(Improve ~ TX + Hospital, weights = Freq, data = d1.rev.df,
                family = binomial(link = "logit"))
d1.revResults = logis.mod.results(d1.modrev)

d2.modrev = glm(Improve ~ TX + Hospital, weights = Freq, data = d2.rev.df,
                family = binomial(link = "logit"))
d2.revResults = logis.mod.results(d2.modrev)

d3.modrev = glm(Improve ~ TX + Hospital, weights = Freq, data = d3.rev.df,
                family = binomial(link = "logit"))
d3.revResults = logis.mod.results(d3.modrev)

d1.revResults
d2.revResults
d3.revResults




# Interaction models
d1.int = glm(Improve ~ TX*Hospital, weights = Freq, data = d1.rev.df,
             family = binomial(link = "logit"))
d1.intResults = logis.mod.results(d1.int)

d2.int = glm(Improve ~ TX*Hospital, weights = Freq, data = d2.rev.df,
             family = binomial(link = "logit"))
d2.intResults = logis.mod.results(d2.int)

d3.int = glm(Improve ~ TX*Hospital, weights = Freq, data = d3.rev.df,
             family = binomial(link = "logit"))
d3.intResults = logis.mod.results(d3.int)

d1.intResults
d2.intResults
d3.intResults




# Pretty num the p-values
pvalPretty = function(pvals){
  prettyVals = prettyNum(pvals,scientific = FALSE)
  return(matrix(prettyVals))
}
pvalPretty(d1.intResults$P.value)
pvalPretty(d2.intResults$P.value)
pvalPretty(d3.intResults$P.value)



# Wald test on models - NOT CORRECT
# anova(d1.modrev, update(d1.modrev, ~1), test = "Chisq")
# anova(d2.modrev, update(d2.modrev, ~1), test = "Chisq")
# anova(d3.modrev, update(d3.modrev, ~1), test = "Chisq")




# Summary of output
d1.cc
d2.cc
d3.cc

BreslowDayTest(d1.cmh.tabs)
BreslowDayTest(d2.cmh.tabs)
BreslowDayTest(d3.cmh.tabs)

mantelhaen.test(d1.cmh.tabs)
mantelhaen.test(d2.cmh.tabs)
mantelhaen.test(d3.cmh.tabs)

d1.revResults
d2.revResults
d3.revResults
pvalPretty(d3.revResults$P.value)

d1.intResults
d2.intResults
d3.intResults
pvalPretty(d3.intResults$P.value)




# library(car) - the p-values are the same as the z p-values already in the model! Yay!
Anova(d1.mod, type = "II", test = "Wald")
Anova(d2.mod, type = "II", test = "Wald")
Anova(d3.mod, type = "II", test = "Wald")





# WEEK 8 HW KEY




###################################################################################
###################################################################################

# Week 9 - More logistic regression

# LOGISTIC REGRESSION AND LOG-BINOMIAL REGRESSION ARE TWO DIFFERENT THINGS!!!
# WHICH METHOD DOES R USE???
# WHAT DOES family = binomial(link = "logit") MEAN??? 
# FROM STACKOVERFLOW IT LOOKS LIKE family = binomial(link = "log") MEANS LOG-BINOMIAL.
## SOURCE: https://stackoverflow.com/questions/25519343/why-are-log-binomial-regression-results-different-in-r-and-sas

# RR vs OR: link = "log" vs link = "logit"




# Seatbelt dataset
car = 
  "1 1 1   34
1 1 0   273
1 0 1   87
1 0 0   221
0 1 1    6
0 1 0   37
0 0 1   12
0 0 0   51"

car = 
  myDATstring.to.df(car,"seat belt injury n")
car
car$seat[car$seat == 1] = "Front"
car$seat[car$seat == 0] = "Back"
car$belt[car$belt == 1] = "With"
car$belt[car$belt == 0] = "Without"
car$injury[car$injury == 1] = "Injured"
car$injury[car$injury == 0] = "Safe"
car

car.tabs = xtabs(n ~ .,
                 data = car)
car.tabs

car.desc = 
  Desc(car.tabs)

car.mar = 
  myDescToAlmightyDF(car.desc,"Injured","Safe","Seat")

# CC car.mar
car.cc = car.mar[,-c(6,13:16)]
car.cc

# Cohort car.mar
car.cohort = car.mar[,-c(7:12)]
car.cohort  




# In class dataset - alcohol and esophageal cancer
esoph = read.sas7bdat("/Users/jamescutler/Desktop/Frequency_Data/wk9_hwsasfile.sas7bdat")

# Don't forget to change factors to characters before creating new dummy variables
esoph$AGE = as.character(esoph$AGE)
esoph$ETOH = as.character(esoph$ETOH)

# Alcohol use dummy variables
esoph$ETOH[esoph$ETOH %in% c("High","Mod")] = "USE"
esoph$ETOH[esoph$ETOH %in% c("Occ","No")] = "NO"

# Age dummy variables - 25-34; 35-54; 55-74; 75+, or ...
# ...                   25-44; 45-54; 55-64; 65+
esoph$AGE[esoph$AGE %in% c("35-44","45-54")] = "35-54"
esoph$AGE[esoph$AGE %in% c("55-64","65-74")] = "55-74"
unique(esoph$AGE)

# Is this dummy grouping of age reasonable?
edt = data.table(esoph)
a = edt[,sum(N),by=AGE]
barplot(a$V1,names.arg = a$AGE,
        main = "Age groups 25-34; 35-54; 55-74; 75+",
        ylab = "Raw counts",xlab = "Age groups")

# Or is the second way better?
esoph = read.sas7bdat("/Users/jamescutler/Desktop/Frequency_Data/wk9_hwsasfile.sas7bdat")
esoph$AGE = as.character(esoph$AGE)
esoph$AGE[esoph$AGE %in% c("25-34","35-44")] = "25-44"
esoph$AGE[esoph$AGE %in% c("65-74","75+")] = "65+"
unique(esoph$AGE)
edt2 = data.table(esoph)
b = edt2[,sum(N),by=AGE]
barplot(b$V1,names.arg = b$AGE,
        main = "Age groups 25-44; 45-54; 55-64; 65+",
        xlab="Age groups",ylab="Raw counts")
# THIS AGE GROUPING LOOKS BETTER - MORE UNIFORM GROUP SIZES

# Now redo the ETOH dummy variable step
esoph$ETOH = as.character(esoph$ETOH)
# Alcohol use dummy variables
esoph$ETOH[esoph$ETOH %in% c("High","Mod")] = "Use_Alcohol"
esoph$ETOH[esoph$ETOH %in% c("Occ","No")] = "Don't_Use"
esoph

edt = data.table(esoph)
new.e = edt[,sum(N),by=c("AGE","ETOH","CASE")]

# Change CASE labels to "Cancer" "No.Cancer"
new.e$CASE[new.e$CASE == 1] = "Cancer"
new.e$CASE[new.e$CASE == 0] = "No.Cancer"

# And make CASE a factor
new.e$CASE = factor(new.e$CASE)


# Logistic regression
mod1 = glm(CASE ~ ETOH + AGE, 
           weights = V1, 
           data = new.e,
           family = binomial(link = "logit"))
summary(mod1)
logis.mod.results(mod1)

# But this is the same problem as before - I still need to relevel the variables to 
## make sure the reference groups are all right. It looks like the only one I need to
## switch here is the outcome, CASE ("Cancer" vs "No.Cancer").
new.e = 
  new.e %>%
  mutate(CASE = relevel(CASE, ref = "No.Cancer"))

# THIS IS THE CORRECT MODEL - THE PREVIOUS ONE WAS JUST TO REMEMBER YOU NEED 
## TO SWITH THE REFERENCE LEVEL IF IT'S NOT ALREADY ALPHABETICALLY IN THE RIGHT SPOT!
mod2 = glm(CASE ~ ETOH + AGE,
           weights = V1,
           data = new.e,
           family = binomial(link = "logit"))

mod2.results = logis.mod.results(mod2)

pvalPretty(mod2.results$P.value)

# Descriptive statistics
# Stacked bar charts with more than 2 confounder categories (facet grid)
## Use dplyr to create percentages by age and case status:
## Source: https://stackoverflow.com/questions/24576515/relative-frequencies-proportions-with-dplyr

# facet = age, x = outcome (cuz CC), fill = exposure
e.props1 = 
  new.e %>%
  group_by(AGE) %>%
  mutate("percent" = V1/sum(V1))

# This is the one (each outcome status within each age group adds it's two percents to 1)
## THIS IS GORGEOUSLY SIMPLE, INTUITIVE AND POWERFUL CODE
e.props2 = 
  new.e %>%
  group_by(AGE,CASE) %>%
  mutate("percent" = V1/sum(V1))

# Make the variable levels more publication-ready
e.props2$CASE = gsub("\\."," ",e.props2$CASE)
e.props2$ETOH = gsub("\\_", " ",e.props2$ETOH)      # Doesn't matter because the legend labels are polished below in scale_fill_manual

# Voila.
ggplot(e.props2, aes(x=CASE,y=percent,fill=ETOH,cumulative=TRUE)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c("#F8766D","#00BFC4"),
                    labels = c("No alcohol","High alcohol")) +
  geom_text(aes(label = paste0(round(percent*100,1),"%")),
            position = position_stack(vjust = .5)) +
  facet_grid(~AGE) +
  labs(title = "Percentage of people who used alcohol, by outcome",
       subtitle = "Stratified by age group",
       x="",y="",fill = "Alcohol consumption")

# How about a plot with age groups by cancer status
new.e = data.table(new.e)
age.e = new.e[,sum(V1),by=c("CASE","AGE")][order(CASE,decreasing = TRUE)]
age.props=
  age.e %>%
  group_by(CASE) %>%
  mutate("percent" = V1/sum(V1))

ggplot(age.props, aes(x=CASE,y=percent,fill=AGE,cumulative=TRUE)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  geom_text(aes(label = paste0(round(percent*100,1),"%")),
            position = position_stack(vjust = .5)) +
  labs(title = "Percentage belonging to each age group, by outcome",
       x="",y="",fill = "Age group")

# And now how about the reverse? Well, because cancer status was fixed, we shouldn't, right?

# Lo and behold the interactions aren't even significant
mod.int = glm(CASE ~ ETOH*AGE,
              weights = V1,
              data = new.e,
              family = binomial(link = "logit"))

int.results = logis.mod.results(mod.int)
pvalPretty(int.results$P.value)


# Final list of analysis outputs
mod2.results
pvalPretty(mod2.results$P.value)
int.results
pvalPretty(int.results$P.value)


new.e
e.props2
age.props
age.tabs =
  xtabs(V1 ~ CASE+AGE,
        data = age.props)

# This doesn't keep the matrix form
(chisq.residuals(age.tabs, std = TRUE))

# This does. You have to convert the table to a matrix before you use tidy to 
## display the residuals. Then tidy will keep it in matrix form.
age.tabs = as.matrix(age.tabs)
tidy(chisq.residuals(age.tabs,std = TRUE))




# A function for converting ctsdf to dtdf (see below for explanation)
countsToCases <- function(df, countcol = "Freq") {
  # Get the row indices to pull from df
  idx <- rep.int(seq_len(nrow(df)), df[[countcol]])
  
  # Drop count column
  df[[countcol]] <- NULL
  
  # Get the rows from df
  df[idx, ]
}
kidsLong = countsToCases(kids,countcol = "count")
# THE THREE MAIN FORMS OF DATA
## 1) Data frame with repeated cells and no counts column: de-tabled df
## 2) Data frame with a counts column: counts df
## 3) Contingency table: table

# So,
## 1 = de-tabled df (dtdf)
## 2 = counts df    (ctsdf)
## 3 = table        (table)

# Conversions from one to another form
## 1 to 2: table() converts dtdf to table; data.frame() converts table to ctsdf
## 1 to 3: table() converts dtdf to table
## 2 to 1: countsToCases() converts ctsdf to dtdf
## 2 to 3: xtabs() converts ctsdf to table
## 3 to 1: data.frame() converts table to ctsdf; countsToCases() converts ctsdf to dtdf
## 3 to 2: data.frame() converts table to ctsdf







###################################################################################
###################################################################################

# Week 10 - more logistic regression
## LR vs Wald vs Score tests
## multicollinearity

# Multicollinearity
## SEE biostats_II script for that


# in class exercise
-3.84 + 1.71
exp(1.71)

a = -3.84
b1 = 1.71
a2 = 1.97
a3 = 2.49
a4 = 2.75
a5 = 2.11

d = exp(a + b1 + a2)
e = exp(a + a2)
d/e
exp(a + b1)/(1+.1188373)




# in class Challenger exercise
chall = 
  "1 66 0  
2 70 1  
3 69 0  
4 68 0  
5 67 0  
6 72 0
7 73 0  
8 70 0  
9 57 1 
10 63 1 
11 70 1 
12 78 0
13 67 0 
14 53 1 
15 67 0 
16 75 0 
17 70 0 
18 81 0
19 76 0 
20 79 0 
21 75 1 
22 76 0 
23 58 1  
24 31 1"
chall = unlist(strsplit(chall,"\n"))
chall = gsub("\\s+$","",chall)
chall = unlist(strsplit(chall," "))

temp = chall[seq(2,length(chall),3)]
distress = chall[seq(3,length(chall),3)]
chall = data.frame(temp = temp,
                   distress = distress,
                   stringsAsFactors = FALSE)

# REMEMBER THIS CODE
chall = data.frame(lapply(chall,as.numeric))

# Since it's univariate, plot the logistic regression model!
ggplot(chall, aes(temp,distress)) +
  geom_point(alpha = .4) +
  stat_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = FALSE)
# This leads me to conclude that as the temperature INCREASES, the probability of 
## distress DECREASES.


#######################################################################################
# NOT SO FAST. THE GLM FUNCTION CAN'T WORK WITH CHARACTER OUTCOMES UNLESS YOU HAVE    #  
## COUNT DATA WITH A weights = n ARGUMENT?                                            #
# Change the enigmatic 0 and 1 to No and Yes now that I've exhausted their usefulness #
## in the plot above.                                                                 #
# chall$distress[chall$distress == 0] = "No"                                          # 
# chall$distress[chall$distress == 1] = "Yes"                                         #
# "No" will be the baseline, or reference group, given that it's alphabetically first.#
## No relevel required. (ALWAYS MAKE THE NON-TERRIBLE, OR, NO-CHANGE, OUTCOME THE     # 
## BASELINE).                                                                         #
#######################################################################################


# Logistic regression model
modchall = glm(distress ~ temp, 
               data = chall,
               family = binomial(link = "logit"))
summary(modchall)
chall.results = logis.mod.results(modchall)
# And sure enough, as my plot indicated, the probability of distress does decrease
## as temperature increases.



# Model diagnostics: Was the model significant?
## Source: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/anova.glm.html

# LR test
LRtest.results = anova(modchall)
LRtest.stat = LRtest.results$Deviance
LRtest.pvalue = 1-pchisq(LRtest.stat,1)
# The Likelihood Ratio test results
LRtest.stat          # X^2 = 10.23671
LRtest.pvalue        # p-value = 0.001376731
LRtest.results$Df    # DF = 1
# BTW ... 
anova(modchall, test = "LRT")
# ... also works--it gives the p-value so you don't have to type another line of code
## to get it.

# The Score test (Lagrange Multiplier test) results
anova(modchall, test = "Rao")
## Score, or Rao, X^2 = 7.9335
## P-value = 0.004853
## DF = 1 of course

# Wald test
## requires: library(car)
Anova(modchall, type = "II", test = "Wald")
# Wald X^2 = 4.6304
# P-value = 0.03141, same as the p-value for temp in the model itself
# DF = 1 of course


# ROC 
# install.packages("pROC")
# requires library(pROC)
challROC = roc(distress ~ temp, data = chall)
plot(challROC, col = "red")
challROC
# Area under the curve: 0.8086


# Hosmer-Lemeshow GOF test
## Null hypothesis: the current model fits well.
# requires library(ResourceSelection)
modchall$y
fitted(modchall)
hoslem.test(modchall$y,fitted(modchall),g=10)
# The resulting p-value, 0.2858, is not quite the same as the one in the key (0.23).
## I have no idea why this is.



# When the Challenger shuttle launched, the temperature was 31 degrees. What was
## the chance of an O-ring distress incident (failure)?
newdata = data.frame(temp = 31)
predict(modchall, newdata = newdata, type = "response")
# The probability of failure was .9996








###################################################################################
###################################################################################

# Week 10 HW
# pros = [big fat txt file dataset]
# pr = myDATstring.to.df(pros,"ID capsule age race dpros dcaps psa vol gleason")
# write.csv(pr, "/Users/jamescutler/Desktop/Frequency_Data/Prostate_Frequency_Class.csv")
pr = read.csv("/Users/jamescutler/Desktop/Frequency_Data/Prostate_Frequency_Class.csv")

class(pr$capsule)
class(pr$dpros)
class(pr$psa)


# Logistic regression - interaction model - DPROS NUMERIC
mod.int = glm(capsule ~ psa*dpros,
              data = pr,
              family = binomial(link = "logit"))
summary(mod.int)
modint.results = logis.mod.results(mod.int)
modint.results
pvalPretty(modint.results$P.value)
# The interaction was only marginally significant (P=.066).


# Logistic regression - main effects model - DPROS NUMERIC
mod.ma = glm(capsule ~ psa + dpros,
             data = pr,
             family = binomial(link = "logit"))
summary(mod.ma)
modma.results = logis.mod.results(mod.ma)
modma.results
pvalPretty(modma.results$P.value)


# Dummy variables for a categorical version of DPROS
pr$dprosCAT = NA
pr$dprosCAT[pr$dpros == 1] = "A_No"
pr$dprosCAT[pr$dpros == 2] = "B_Left"
pr$dprosCAT[pr$dpros == 3] = "C_Right"
pr$dprosCAT[pr$dpros == 4] = "D_Bilobar"
# Logistic regression - interaction model - DPROS CATEGORICAL
modint.cat = glm(capsule ~ psa*dprosCAT,
                 data = pr,
                 family = binomial(link = "logit"))
summary(modint.cat)
modintcat.results = logis.mod.results(modint.cat)
modintcat.results
# Interaction was not significant.


# Logistic regression - main effects - DPROS CATEGORICAL
modma.cat = glm(capsule ~ psa + dprosCAT,
                data = pr,
                family = binomial(link = "logit"))
summary(modma.cat)
modmacat.results = logis.mod.results(modma.cat)
modmacat.results
pvalPretty(modmacat.results$P.value)
# All variables and factor levels were significant.


# Diagnostics
# LRT
# anova(modma.cat, test = "LRT")                # Doesn't give an "overall p-value for the model"
## Source: https://rcompanion.org/rcompanion/e_06.html - VERY HELPFUL!
anova(modma.cat, update(modma.cat, ~1), test = "LRT")   # chisq = 82.14,  df = 4
1-pchisq(82.14,4)
lmtest::lrtest(modma.cat)                               # chisq = 82.14,  df = 4

# Score
# anova(modma.cat, test = "Rao")                # Doesn't give an "overall p-value for the model"
anova(modma.cat, update(modma.cat, ~1), test = "Rao")   # chisq = 69.596, df = 4
1-pchisq(69.596,4)

# Wald
# Anova(modma.cat, type = "II", test = "Wald")  # Doesn't give an "overall p-value for the model" 
# anova(modma.cat, update(modma.cat, ~1), test="Chisq") # GIVES SAME RESULT AS LRT
# Anova(modma.cat, update(modma.cat, ~1), type = "II", test = "Wald") # DOESN'T WORK
# lmtest::waldtest(modma.cat)                             # F = 13.123,     df = 4
# 1-pchisq(13.123,4)                                      # So what is it?  
# 1-pf(13.123,df1=4,df2=4)                                # So what is it?
lmtest::waldtest(modma.cat, test = "Chisq")             # chisq = 52.49,  df = 4
1-pchisq(52.494,4)



# ROC
library(pROC)
prosROC = roc(capsule ~ psa + dpros, data = pr)
prosROC
plot(prosROC$psa, col = "red")
plot(prosROC$dpros, col = "green", add = TRUE)
auc(prosROC$psa)
auc(prosROC$dpros)


# Hosmer-Lemeshow
## Null: the current model fits well.
ResourceSelection::hoslem.test(modma.cat$y, fitted(modma.cat),g=10)



# Write the model in terms of the coefficients (see HW assignment)
summary(modma.cat)

# Report the ORs, CIs, P-values
modmacat.results
pvalPretty(modmacat.results$P.value)

# Report the ROC and Hoslem-Lemeshow results
auc(prosROC$psa)
auc(prosROC$dpros)

# Final statistical and practical interpretations
summary(modma.cat)
modmacat.results




# Make a prediction - THE FUN PART!!!

## What should my new PSA value be?
hist(pr$psa, breaks = 30)
ggplot(pr, aes(psa,capsule)) +
  geom_point(alpha = .4) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = FALSE)
# How about 60 mg/ml

## What should my new DPROS subcategory be?
hist(pr$dpros)
# How about left side

# If a patient walks in and gets examined, and they have a PSA level of 60 mg/ml,
## with a nodule on the left side, what are their chances of capsule penetration?
newdata = data.frame(psa = 60, dprosCAT = "B_Left")
predict(modma.cat, newdata = newdata, type = "response")
# Answer: They have an 86% chance of capsule penetration.

# PSA of 60, with right side nodule
newdata = data.frame(psa = 60, dprosCAT = "C_Right")
predict(modma.cat, newdata = newdata, type = "response")
# Answer: 92.8% chance of penetration.

# PSA of 60, with bilobar nodule
newdata = data.frame(psa = 60, dprosCAT = "D_Bilobar")
predict(modma.cat, newdata = newdata, type = "response")
# Answer: Still 92.8% chance of penetration.


# PSA of 40, with left side nodule
newdata = data.frame(psa = 40, dprosCAT = "B_Left")
predict(modma.cat, newdata = newdata, type = "response")
# Answer: 69.3% chance of penetration.

# PSA of 40, with right side nodule
newdata = data.frame(psa = 40, dprosCAT = "C_Right")
predict(modma.cat, newdata = newdata, type = "response")
# Answer: 82.5% chance of penetration.


# PSA of 20, with left side nodule
newdata = data.frame(psa = 20, dprosCAT = "B_Left")
predict(modma.cat, newdata = newdata, type = "response")
# Answer: 45.3% chance of penetration.

# PSA of 10, with left side nodule
newdata = data.frame(psa = 10, dprosCAT = "B_Left")
predict(modma.cat, newdata = newdata, type = "response")
# Answer: 33.4% chance of penetration.


# PSA of 20, with no nodule
newdata = data.frame(psa = 20, dprosCAT = "A_No")
predict(modma.cat, newdata = newdata, type = "response")
# Answer: 23.3% chance of penetration.

# PSA of 10, with no nodule
newdata = data.frame(psa = 10, dprosCAT = "A_No")
predict(modma.cat, newdata = newdata, type = "response")
# Answer: 15.6% chance of penetration.

# PSA of 10, with right side nodule
newdata = data.frame(psa = 10, dprosCAT = "C_Right")
predict(modma.cat, newdata = newdata, type = "response")
# Answer: 51.2% chance of penetration.









###################################################################################
###################################################################################

# Week 11 in class exercise

kids = 
  "1	1	1	1	4
1	1	1	0	349
1	1	2	1	2 
1	1	2	0	232
1	1	3	1	8
1	1	3	0	166
1	1	4	1	4
1	1	4	0	48
1	2	1	1	13
1	2	1	0	64
1	2	2	1	27
1	2	2	0	84
1	2	3	1	47
1	2	3	0	91
1	2	4	1	39
1	2	4	0	57
2	1	1	1	9
2	1	1	0	207
2	1	2	1	7
2	1	2	0	201
2	1	3	1	6
2	1	3	0	120
2	1	4	1	5
2	1	4	0	47
2	2	1	1	33
2	2	1	0	72
2	2	2	1	64
2	2	2	0	95
2	2	3	1	74
2	2	3	0	110
2	2	4	1	123
2	2	4	0	90
3	1	1	1	12
3	1	1	0	126
3	1	2	1	12
3	1	2	0	115
3	1	3	1	17
3	1	3	0	92
3	1	4	1	9
3	1	4	0	41
3	2	1	1	38
3	2	1	0	54
3	2	2	1	93
3	2	2	0	92
3	2	3	1	148
3	2	3	0	100
3	2	4	1	224
3	2	4	0	65
4	1	1	1	10
4	1	1	0	67
4	1	2	1	17
4	1	2	0	79
4	1	3	1	6
4	1	3	0	42
4	1	4	1	8
4	1	4	0	17
4	2	1	1	49
4	2	1	0	43
4	2	2	1	119
4	2	2	0	59
4	2	3	1	198
4	2	3	0	73
4	2	4	1	414
4	2	4	0	54"

kids = myDATstring.to.df(kids,"iq parent ses college count")
kids



mod1.int = glm(college ~ iq*parent + ses, weights = count,
               data = kids, family = binomial(link = "logit"))
mod1.int.results = logis.mod.results(mod1.int)
pvalPretty(mod1.int.results$P.value)

mod2.int = glm(college ~ iq*ses + parent, weights = count,
               data = kids, family = binomial(link = "logit"))
mod2.int.results = logis.mod.results(mod2.int)
mod2.int.results

mod3.int = glm(college ~ iq + parent*ses, weights = count,
               data = kids, family = binomial(link = "logit"))
mod3.int.results = logis.mod.results(mod3.int)
mod3.int.results

combn(c("iq","parent","ses"),2)

mod.me = glm(college ~ iq + parent + ses, weights = count,
             data = kids, family = binomial(link = "logit"))
mod.me.results = logis.mod.results(mod.me)
mod.me.results
kids$catses = NA
kids$catses[kids$ses == 1] = "A"
kids$catses[kids$ses == 2] = "B"
kids$catses[kids$ses == 3] = "C"
kids$catses[kids$ses == 4] = "D"

kids$catiq = NA
kids$catiq[kids$iq == 1] = "A"
kids$catiq[kids$iq == 2] = "B"
kids$catiq[kids$iq == 3] = "C"
kids$catiq[kids$iq == 4] = "D"
mod.me2 = glm(college ~ catiq + parent + catses, weights = count,
              data = kids, family = binomial(link = "logit"))
modme2.results = logis.mod.results(mod.me2)
modme2.results
class(kids$ses)
unique(kids$parent)


# LRT
anova(mod.me2, update(mod.me2, ~1), test = "LRT")

# Score
anova(mod.me2, update(mod.me2, ~1), test = "Rao") 

# Wald
lmtest::waldtest(mod.me2, test = "Chisq")    


# write the model
mod.me2


# p-values, ORs, CIs
modme2.results
pvalPretty(modme2.results$P.value)


# ROC
countsToCases <- function(df, countcol = "Freq") {
  # Get the row indices to pull from df
  idx <- rep.int(seq_len(nrow(df)), df[[countcol]])
  
  # Drop count column
  df[[countcol]] <- NULL
  
  # Get the rows from df
  df[idx, ]
}
kidsLong = countsToCases(kids,countcol = "count")
kidsROC = roc(college ~ iq + parent + ses, data = kidsLong)
kidsROC
plot(kidsROC$iq, col = "red")
plot(kidsROC$parent, col = "green", add = TRUE)
plot(kidsROC$ses, col = "blue",add=TRUE)


# Hosmer Lemeshow - NOT NEEDED!






###################################################################################
###################################################################################

# HW 11

# How to report everything:
# 1) gender proportions among treatment groups
# 2) age proportions among treatment groups
# 3) descriptive stats: gender by outcome
# 4) descriptive stats: age by outcome
# 5) descriptive stats: treatment by outcome
# 5) chi-square univariate results
# 6) multivariable interaction models - testing for interaction
# 7) final multivariable model results
# 8) final model dx

fib = read_xls("/Users/jamescutler/Desktop/Frequency_Data/Fibro.xls")
fib

# 1='Female'
# 0='Male';
# value Tx
# 1='Active Treatment'
# 0='Placebo';
# value age
# 1='Age Below 50'
# 0='Greater or equal 50';
# value decrease
# 1='Yes'
# 0='No';



fibLong = countsToCases(fib,countcol = "count")
fibLong = as.data.frame(fibLong)
summary(fibLong)


#################### FOR SUMMARY STATS #####################
fibLong2 = fibLong
fibLong2 = apply_labels(fibLong,                           #
                        gender = "Gender",                 #
                        tx = "Treatment",                  #
                        age50 = "Age group, yrs",          #
                        decrease = "Pain decrease (Y/N)")  #
fibLong2$gender[fibLong2$gender == 0] = "Male"             #
fibLong2$gender[fibLong2$gender == 1] = "Female"           #

fibLong2$tx[fibLong2$tx == 0] = "Placebo"                  #
fibLong2$tx[fibLong2$tx == 1] = "Active"                   #

fibLong2$age50[fibLong2$age50 == 0] = "50 or greater"      #
fibLong2$age50[fibLong2$age50 == 1] = "Under 50"           #

fibLong2$decrease[fibLong2$decrease == 0] = "No"           #
fibLong2$decrease[fibLong2$decrease == 1] = "Yes"          #
############################################################


# Univariate - NOT JUST 3 CHI SQUARE TESTS--5! 
# YOU DON'T JUST TEST FOR ASSOCIATION BETWEEN EACH COVARIATE AND THE OUTCOME,
# BUT ALSO BETWEEN EACH DEMOGRAPHIC COVARIATE AND TX!!
combn(c("gender","tx","age","outcome"),2)
## gender and outcome
gender.tab = table(fibLong2$gender,fibLong2$decrease)
gender.tab
chisq.test(gender.tab)
## tx and outcome
tx.tab = table(fibLong2$tx,fibLong2$decrease)
tx.tab
chisq.test(tx.tab)
## age and outcome
age.tab = table(fibLong2$age50,fibLong2$decrease)
age.tab
chisq.test(age.tab)

## age and tx
chisq.test(table(fibLong2$age50,fibLong2$tx),correct = FALSE)
## gender and tx
chisq.test(table(fibLong2$gender,fibLong2$tx),correct = FALSE)



# Main effects
modME = glm(decrease ~ age50 + tx + gender, weights = count, data = fib,
            family = binomial(link = "logit"))
modmeR = logis.mod.results(modME)
modmeR


# Interaction
## All 3 interactions
modallint = glm(decrease ~ age50*tx*gender, weights = count, data = fib,
                family = binomial(link = "logit"))
modallintR = logis.mod.results(modallint)
modallintR
pvals = pvalPretty(modallintR$P.value)
modallintR$P.value = pvals
modallintR
# AGE AND GENDER HAVE A SIGNIFICANT INTERACTION

## tx and age
modTA = glm(decrease ~ age50*tx + gender, weights = count, data = fib,
            family = binomial(link = "logit"))
modtaR = logis.mod.results(modTA)
modtaR
pvalPretty(modtaR$P.value)
## tx and gender
modTG = glm(decrease ~ age50 + tx*gender, weights = count, data = fib,
            family = binomial(link = "logit"))
modtgR = logis.mod.results(modTG)
modtgR
pvalPretty(modtgR$P.value)
## age and gender - STILL SIGNIFICANT
modAG = glm(decrease ~ age50*gender + tx, data = fibLong,
            family = binomial(link="logit"))
modagR = logis.mod.results(modAG)
modagR
pvals = pvalPretty(modagR$P.value)
modagR$P.value = pvals
modagR

# Since there is an interaction between age and gender, we could stratify the 
## models by gender, and then also by age - AMAZING!!!!!!!!
## Source: https://rpubs.com/corey_sparks/335463
modMale = glm(decrease ~ age50 + tx, data = fibLong,
              family = binomial(link="logit"), subset = gender==0)
modMaleR = logis.mod.results(modMale)
modMaleR # Among males, the odds of a pain decrease are 50% lower for the young compared to the old. 

modFem = glm(decrease ~ age50 + tx, data = fibLong,
             family = binomial(link="logit"), subset = gender==1)
modFemR = logis.mod.results(modFem)
modFemR # Among females, the odds of a pain decrease are 7.9 times greater for the young compared to the old. 

# Age
modYoung = glm(decrease ~ gender + tx, data = fibLong,
               family = binomial(link="logit"), subset = age50==1)
modYoungR = logis.mod.results(modYoung)
modYoungR # Among younger patients, the odds of a pain decrease are 9.2 times higher for females compared to males.

modOld = glm(decrease ~ gender + tx, data = fibLong,
             family = binomial(link="logit"), subset = age50==0)
modOldR = logis.mod.results(modOld)
modOldR # Among older patients, the odds of a pain decrease are 48.5% lower for females compared to males.

# But if you want to use LogisticDx::gof on these stratified models, you need to 
## do all of this crap instead! Not only do you have to create two new dataframes,
## you have to change the groups (g=11) in the gof function to something that will
## work with the Hosmer-Lemeshow part of the source code.
fibMale = fibLong[fibLong$gender == 0,]
fibFemale = fibLong[fibLong$gender == 1,]
modMale2 = glm(decrease ~ age50 + tx, data = fibMale,
               family = binomial(link="logit"))
modMale2R = logis.mod.results(modMale2)
modFem2 = glm(decrease ~ age50 + tx, data = fibFemale,
              family = binomial(link="logit"))
gof(modMale2,g=11)
gof(modFem2,g=11)





# Pearson gof, Deviance, and C statistic all in one! AMAZING!!!!!!!!!!!!
# library(LogisticDx)

# GET A MODEL FROM THE LONG DATAFRAME - THE COUNTS DATAFRAME IS CAUSING PROBLEMS WITH A
## LOT OF THE DIAGNOSTICS PACKAGES
LogisticDx::gof(modAG)       # USE modAG2!!! - DOESN'T WORK WITH modAG APPARENTLY BECAUSE modAG IS FROM THE COUNTS DATAFRAME




# Now for summarizing descriptive statistics and the gender-age interaction
fibdt = data.table(fibLong2)
fibdt[,.N,by="gender"]
fibdt
d = fibdt[,.N,by=c("tx","gender")][order(tx)]
fibdt[,.N,by="age50"]

barplot(d$N, names.arg = d$gender)
table(fibLong2$gender,fibLong2$tx)
barplot(table(fibLong2$gender,fibLong2$tx))

# Do the treatment groups have equal gender proportions?
prop.test(table(fibLong2$gender,fibLong2$tx))
prop.test(table(fibLong2$age50,fibLong2$tx))
barplot(table(fibLong2$age50,fibLong2$tx))

gtx = table(fibLong2$gender,fibLong2$tx)
chisq.residuals(gtx, std = TRUE)

gpain = table(fibLong2$gender,fibLong2$decrease)

prop.table(gpain,1)


# Percentages by treatment group
g.nums = fibdt[,.N,by=c("tx","gender","decrease")][order(tx)]

g.props = 
  g.nums %>%
  group_by(tx,gender) %>%
  mutate("percent" = N/sum(N))
g.props


# Proportions of each gender across age groups
ga.nums = fibdt[,.N,by=c("age50","gender")][order(gender)]
ga.nums  
table(fibLong2$age50)

ga.props = 
  ga.nums %>%
  group_by(gender) %>%
  mutate("percent" = N/sum(N))
ga.props





# The OR table

# ORs
## tx
## age - female
## age - male
## gender - young
## gender - old
ORs = 
  rbind(
    modagR[3,1:3],
    modFemR[1,1:3],
    modMaleR[1,1:3],
    modYoungR[1,1:3],
    modOldR[1,1:3]
  )
rownames(ORs) = c("Treatment (active vs placebo)",
                  "Female age (young vs old)",
                  "Male age (young vs old)",
                  "Young gender (female vs male)",
                  "Old gender (female vs male)")
save(ORs, file = "/Users/jamescutler/Desktop/Frequency_Data/Report3_ORs_Table")
# Use this in RMarkdown to knit the HTML table

# Models
#modAG
modagR

#modMale
malePvals = pvalPretty(modMaleR$P.value)
modMaleR$P.value = malePvals
modMaleR

#modFem
modFemR

#modYoung
modYoungR

#modOld
modOldR
pvalPretty(modOldR$P.value)








#######################################################################################
#######################################################################################

# Week 12 - McNemar, Cohen's Kappa, Bowker's symmetry, exact McNemar's

control = c("exposed","unexposed")
case = c("exposed","unexposed")
cc = matrix(c(45,20,
              10,50),2,2,byrow = TRUE)
dimnames(cc) = list("Case" = case,"Control" = control)
cc

# McNemar's OR
20/10

# 95% CI for log(OR)
logup = log(2) + 1.96*sqrt(1/20 + 1/10)
loglow = log(2) - 1.96*sqrt(1/20 + 1/10)

# 95% CI for OR
exp(logup)
exp(loglow)

mcnemar.test(cc, correct = FALSE) # huge difference between continuity correction and no correction


# Load library irr
library(irr)

pathA = c("Stage_0","State_1-3","Stage_4")
pathB = c("Stage_0","State_1-3","Stage_4")
lung = matrix(c(78,5,0,
                6,56,13,
                0,10,32),3,3,byrow = TRUE)
dimnames(lung) = list("Pathologist 1" = pathA, "Pathologist 2" = pathB)
lung

# kappa2(lung) # THINKS IT HAS MORE THAN TWO RATERS

data.frame(lung)


data("anxiety")
head(anxiety)

## 1 = de-tabled df (dtdf)
## 2 = counts df    (ctsdf)
## 3 = table        (table)

# Conversions from one to another form
## 1 to 2: table() converts dtdf to table; data.frame() converts table to ctsdf
## 1 to 3: table() converts dtdf to table
## 2 to 1: countsToCases() converts ctsdf to dtdf
## 2 to 3: xtabs() converts ctsdf to table
## 3 to 1: HECK YEAH: use my tabletoDFforKappa() function - it converts table to dtdf
## 3 to 2: data.frame() converts table to ctsdf

a2 = anxiety[,1:2]
a2
table(a2)

as.data.frame.matrix(lung)

lungdf = data.frame(path1 = c(rep("Stage 0",78+5),
                              rep("Stage 1-3",6+56+13),
                              rep("Stage 4",10+32)),
                    path2 = c(rep("Stage 0",78),
                              rep("Stage 1-3",5),
                              rep("Stage 0",6),
                              rep("Stage 1-3",56),
                              rep("Stage 4",13),
                              rep("Stage 1-3",10),
                              rep("Stage 4",32))
)
lungdf
lungdf$path1 = as.character(lungdf$path1)
lungdf$path2 = as.character(lungdf$path2)

kappa2(lungdf) # works!

lvec = as.vector(t(lung))
lvec


lungdf2 = tabletoDFforKappa(lung,path1,path2)
lungdf2
kappa2(lungdf2)





#######################################################################################
#######################################################################################

# Week 12 HW

path1 = c(paste0("B",1:5))
path2 = c(paste0("B",1:5))
br = matrix(c(7,33 ,2 ,0,0  ,
              2,384,18,0,1  ,
              1,29 ,46,0,14 ,
              0,0  ,0 ,0,1  ,
              0,0  ,0 ,2,225),5,5,byrow = TRUE)
dimnames(br) = list("Pathologist_1" = path1, "Pathologist_2" = path2)
table.margins(br)

# Table to dtdf!
h = tbToDF(br,path1,path2)

kappa2(h)

cl = tbToDF(lung,pathA,pathB)
kappa2(cl)
