---
title: 'carbohydrate, glucose, lactose feeding - peroxide survival - Elises experiemnt from Oct 2017'
author: "Ben"
date: "10/25/2017"
---
  
  # survival analysis of data from a peroxide trial done on four DGRP lines, each fed a range of carbohydrate concentrations (if the flies were avaialble) in their food prior to the peroxide trial:
library(data.table)
library(splitstackshape)
library(dplyr)
library("tidyr")
library('survival')
install.packages('Hmisc') # allows errbar()
library(Hmisc)
library(psych)

mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

setwd("~/Google_Drive/Documents/Peroxide/Elises.Feeding.Experiments")
dat <- read.csv("maltose.glucose.lactose.feeding.expt.csv", header=T, stringsAsFactors=T)
head(dat)
table(dat$UniqueName)
table(dat$IntDeaths, dat$Censored)

# with carried data, I should use IntDeaths as 'Deaths'.
# to make a proper right-censored Surv object, I should dat to make a column of 'event', with 1 for each death
# and 0 for each censored fly.   There should be one row per event and only 0=censored or 1= death (no other entries):

# dLife makes IntDeaths with >1 death per census if >1 death occured, so these need to be replicated and replaced with 1s.  Here's how:
### get rows with >1 death
dat2 <- dat[dat$IntDeaths>1, ]
head(dat2)

### expand to give one row per death per time
dat3 <- expandRows(dat2, "IntDeaths")
head(dat3)
#### Note, the splitstackshape way ends up removing the IntDeaths column.  Need to add back but now all values =1
dat3$IntDeaths <- 1
head(dat3)
head(dat)
dat3 <- dat3[, c(1:9,13,10:12)]
head(dat3)

### add this to the rows with <=1 death
dim(dat[!dat$IntDeaths>1, ])
dat4 <- rbind(dat3, dat[!dat$IntDeaths>1, ])
head(dat4)

table(dat4$IntDeaths, dat4$Censored) # this will check if there are censored flies among the deaths
# eliminate the censored flies in the deaths data:
dat4$Censored <- 0
table(dat4$IntDeaths, dat4$Censored) 
## keep only the rows with a death:
dat4 <- dat4[dat4$IntDeaths == 1, ]
table(dat4$IntDeaths)

### do this same procedure with censored flies######:
### get rows with >1 Censored
dat5 <- dat[dat$Censored>1, ]
head(dat5)
### expand to give one row per censored fly per time
dat6 <- expandRows(dat5, "Censored")
head(dat6)
####### Add back the Censored column, but now all values = 1
dat6$Censored <- 1
head(dat6)
head(dat)
dat6 <- dat6[, c(1:8,13,9:12)]
head(dat6)

### add this to the rows with <=1 Censored
dim(dat4[!dat4$Censored>1, ])
dat7 <- rbind(dat6, dat4[!dat4$Censored>1, ])
head(dat7)

table(dat7$Censored)
# eliminate the deaths from these censored-only data:
dat7$IntDeaths <- 0
table(dat7$Censored, dat7$IntDeaths)
# eliminate rows where censored is not 1
dat7[dat7$Censored == 1, ] -> dat7

# you now have two data frames, one with all of the censored flies, and the other with all of the deaths.  You will rbind these two dfs and this will give a column (IntDeaths) that contains either 1 or 0.  The 1 is when a death occured and 0 is when a censored fly left the experiment.
dat8 <- rbind(dat4, dat7)

#### make a new df$event column ...which should just be the same as IntDeaths because when it is 0 is when a fly was censored.
dat8$event <- dat8$IntDeaths
table(dat8$event)
head(dat8)
dat8 <- separate(dat8, UniqueName, into=c('line', 'carbohydrate', 'treatment'), sep = '/', remove=F)
head(dat8)
str(dat8)
dat8$treatment <- as.factor(dat8$treatment)
dat8[, 'carbohydrate'] <- as.factor(dat8[, 'carbohydrate'])
levels(dat8$carbohydrate)
dat8$carbohydrate <- factor(dat8$carbohydrate, levels = c("water", "maltose", "dextrose", "lactose"))
head(dat8)
str(dat8)
table(dat8$event, dat8$UniqueName) # to get sample sizes

### get survival perameters (e.g. mean, median, etc) with survfit
dat.stats <- print(survfit(Surv(AgeH, event) ~ dat8$UniqueName, data=dat8), print.rmean=T) 
as.data.frame(summary(dat.stats)$table) -> dat.stats
head(dat.stats)
dat.stats$ID <- rownames(dat.stats)
dat.stats <- separate(dat.stats, ID, into=c('line', 'carbohydrate', 'treatment'), sep = '/', remove=F)
head(dat.stats)
dat.stats <- separate(dat.stats, line, into=c('junk', 'line'), sep = '=', remove=T)
head(dat.stats)
dat.stats <- dat.stats[ ,-c(10:11)]
head(dat.stats)
dat.stats <- dat.stats[ ,c(10:12, 1:9)]
head(dat.stats)
rownames(dat.stats)=NULL
colnames(dat.stats) <- c("line", "carbohydrate", 'treatment', "records", "n.max",  "n.start", "event", "rmean", "se.rmean", "median", "0.95LCL", "0.95UCL") ## this step is done to get rid of the astrisks in two of the column names
head(dat.stats)
str(dat.stats)
dat.stats[, 'carbohydrate'] <- as.factor(dat.stats[, 'carbohydrate'])
dat.stats$carbohydrate <- factor(dat.stats$carbohydrate, levels = c("water", "maltose", "dextrose", "lactose"))
dat.stats[, 'treatment'] <- as.factor(dat.stats[, 'treatment'])
str(dat.stats)
table(dat.stats$line)
table(dat.stats$event, dat.stats$treatment)
par(mfrow=c(1,1))
subset(dat.stats, treatment == 'peroxide') -> treats
str(treats)

par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
stripchart(rmean ~ carbohydrate, data=treats[treats$line == 'RAL_161', ], vertical=T, pch=19, col=mycol4[1], cex=2)
stripchart(rmean ~ carbohydrate, data=treats[treats$line == 'RAL_239', ], vertical=T, pch=19, col=mycol4[2], cex=2)
stripchart(rmean ~ carbohydrate, data=treats[treats$line == 'RAL_320', ], vertical=T, pch=19, col=mycol4[3], cex=2)
stripchart(rmean ~ carbohydrate, data=treats[treats$line == 'RAL_627', ], vertical=T, pch=19, col=mycol4[4], cex=2)

### Figure for paper #####
par(mfrow=c(2, 2))
par(mar=c(3.5, 3, 1, 2))

#RAL_239
treats$carbohydrate[treats$line == 'RAL_239'] -> carbs
carbs[c(4, 3, 1, 2)] -> carbs
treats$rmean[treats$line == 'RAL_239'] -> mean
mean[c(4, 3, 1, 2)] -> mean
treats$se.rmean[treats$line == 'RAL_239'] -> se
se[c(4, 3, 1, 2)] -> se

plot(x=1:4, y=mean, col=0, ylim=c(70, 100), xlim=c(.5, 4.5), xaxt='n', ann=F, las=1)
arrows(x0 = 1:4, y0 = mean, x1 = 1:4, y1 = mean + se, angle = 90, length = 0.05, lwd = 1)
arrows(x0 = 1:4, y0 = mean, x1 = 1:4, y1 = mean - se, angle = 90, length = 0.05, lwd = 1)
points(x= 1:4, y = mean, pch = c(19), cex = 2, col=mycol4[2])
text(1:4, par('usr')[3], labels = labels, srt = 45, adj = c(1.1,1.2), xpd = TRUE, cex=.9)

#RAL_320
treats$carbohydrate[treats$line == 'RAL_320'] -> carbs
carbs[c(4, 3, 1, 2)] -> carbs
treats$rmean[treats$line == 'RAL_320'] -> mean
mean[c(4, 3, 1, 2)] -> mean
treats$se.rmean[treats$line == 'RAL_320'] -> se
se[c(4, 3, 1, 2)] -> se

plot(x=1:4, y=mean, col=0, ylim=c(130, 160), xlim=c(.5, 4.5), xaxt='n', ann=F, las=1)
arrows(x0 = 1:4, y0 = mean, x1 = 1:4, y1 = mean + se, angle = 90, length = 0.05, lwd = 1)
arrows(x0 = 1:4, y0 = mean, x1 = 1:4, y1 = mean - se, angle = 90, length = 0.05, lwd = 1)
points(x= 1:4, y = mean, pch = c(19), cex = 2, col=mycol4[3])
text(1:4, par('usr')[3], labels = labels, srt = 45, adj = c(1.1,1.2), xpd = TRUE, cex=.9)

# RAL_627
treats$carbohydrate[treats$line == 'RAL_627'] -> carbs
carbs[c(4, 3, 1, 2)] -> carbs
treats$rmean[treats$line == 'RAL_627'] -> mean
mean[c(4, 3, 1, 2)] -> mean
treats$se.rmean[treats$line == 'RAL_627'] -> se
se[c(4, 3, 1, 2)] -> se

plot(x=1:4, y=mean, col=0, ylim=c(50, 75), xlim=c(.5, 4.5), xaxt='n', ann=F, las=1)
arrows(x0 = 1:4, y0 = mean, x1 = 1:4, y1 = mean + se, angle = 90, length = 0.05, lwd = 1)
arrows(x0 = 1:4, y0 = mean, x1 = 1:4, y1 = mean - se, angle = 90, length = 0.05, lwd = 1)
points(x= 1:4, y = mean, pch = c(19), cex = 2, col=mycol4[4])
text(1:4, par('usr')[3], labels = labels, srt = 45, adj = c(1.1,1.2), xpd = TRUE, cex=.9)

#RAL_161  - omit due to lack of data on glucose condition (insufficient flies to test this condition)
levels(treats$carbohydrate[treats$line == 'RAL_161']) -> carbs
treats$rmean[treats$line == 'RAL_161'] -> mean
mean[3] -> mean[4]
mean[3] <- NA
treats$se.rmean[treats$line == 'RAL_627'] -> se
se[c(4, 3, 1, 2)] -> se
se[3] -> se[4]
se[3] <- NA

plot(x=1:4, y=mean, col=0, ylim=c(65, 105), xlim=c(.5, 4.5), xaxt='n', ann=F)
arrows(x0 = 1:4, y0 = mean, x1 = 1:4, y1 = mean + se, angle = 90, length = 0.05, lwd = 1)
arrows(x0 = 1:4, y0 = mean, x1 = 1:4, y1 = mean - se, angle = 90, length = 0.05, lwd = 1)
points(x= 1:4, y = mean, pch = c(19, 18, 17, 15), cex = 1.5, col=mycol4[1])


##############

par(mfrow=c(1,1))
par(mar=c(6,6,6,6))
plot(rmean ~ numeric.carbohydrate, data=treats[treats$line == 'RAL_161', ], ylim=c(40, 145), col='white', ylab = 'mean lifespan (h)', xlab='supplemental carbohydrate (%)', cex.lab = 2, cex.axis = 1.5)
lines(rmean ~ numeric.carbohydrate, data=treats[treats$line == 'RAL_161',] , type='b', pch=19, lwd=2, col = mycol4[1])
lines(rmean ~ numeric.carbohydrate, data=treats[treats$line == 'RAL_239',] , type='b', pch=19, lwd=2, col = mycol4[2])
lines(rmean ~ numeric.carbohydrate, data=treats[treats$line == 'RAL_320',] , type='b', pch=19, lwd=2, col = mycol4[3])
lines(rmean ~ numeric.carbohydrate, data=treats[treats$line == 'RAL_627',] , type='b', pch=19, lwd=2, col = mycol4[4])

par(xpd=TRUE)
legend(2.2,105 ,c(' RAL 161', ' RAL 239', ' RAL 320', ' RAL 627'), pch = c(19), col=mycol4, bty='n', y.intersp = 1.2)


## code to get error bars from another analysis:
plot(fem.mean.se$genotype, fem.mean.se$line.rmean, pch=19, main='Females', xlab='genotype', ylab='manganese lifespan (h +/- se)', cex.lab=2, cex.main=3, cex.axis=1.5)
errbar(fem.mean.se$genotype, fem.mean.se$line.rmean, fem.mean.se$line.rmean + fem.mean.se$line.se, fem.mean.se$line.rmean - fem.mean.se$line.se, add=T, cap=0.003)
############

### test for the effect of different levels of carbohydrate:
aov(rmean ~ carbohydrate + line, data=treats) -> aov.dat.line
summary(aov.dat.line)
TukeyHSD(x=aov.dat.line, 'carbohydrate')  #get pairwise tests between the carbohydrate levels along with adjusted pValues
par(mfrow=c(1,1))
plot(aov.dat.line) # hit enter in the console below to see some diagnostic plots for the linear model

# use log-rank test to test for effec of maltose, glucose or lactose on lifespan compared to control:
head(dat8)
subset(dat8, treatment=='peroxide') -> temp
head(temp)

i = 'RAL_627'
subset(temp, line == i) -> terp
print(table(terp$carbohydrate))
survdiff(Surv(AgeH, event) ~ UniqueName, data= terp[ terp$carbohydrate == 'water' |  terp$carbohydrate == 'lactose', ])

i='RAL_320'
subset(temp, line == i) -> terp
print(table(terp$carbohydrate))
survdiff(Surv(AgeH, event) ~ UniqueName, data= terp[ terp$carbohydrate == 'water' |  terp$carbohydrate == 'lactose', ])

i='RAL_239'
subset(temp, line == i) -> terp
print(table(terp$carbohydrate))
survdiff(Surv(AgeH, event) ~ UniqueName, data= terp[ terp$carbohydrate == 'water' |  terp$carbohydrate == 'lactose', ])

