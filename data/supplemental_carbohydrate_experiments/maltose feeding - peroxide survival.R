---
title: 'maltose feeding - peroxide survival'
author: "Ben"
date: "10/4/2017"
output: html_document
---

# survival analysis of data from a peroxide trial done on four DGRP lines, each fed a range of maltose concentrations (if the flies were avaialble) in their food prior to the peroxide trial:
library(data.table)
library(splitstackshape)
library(dplyr)
library("tidyr")
library('survival')
install.packages('Hmisc') # allows errbar()
library(Hmisc)

mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")

malt <- read.csv("~/Google_Drive/Documents/Peroxide/Elises.Feeding.Experiments/maltose.trial.1.csv", header=T, stringsAsFactors=T)
head(malt)
table(malt$UniqueName)
table(malt$IntDeaths, malt$Censored)

# with carried data, I should use IntDeaths as 'Deaths'.
# to make a proper right-censored Surv object, I should malt to make a column of 'event', with 1 for each death
# and 0 for each censored fly.   There should be one row per event and only 0=censored or 1= death (no other entries):

# dLife makes IntDeaths with >1 death per census if >1 death occured, so these need to be replicated and replaced with 1s.  Here's how:
### get rows with >1 death
malt2 <- malt[malt$IntDeaths>1, ]
head(malt2)

### expand to give one row per death per time
malt3 <- expandRows(malt2, "IntDeaths")
head(malt3)
#### Note, the splitstackshape way ends up removing the IntDeaths column.  Need to add back but now all values =1
malt3$IntDeaths <- 1
head(malt3)
head(malt)
malt3 <- malt3[, c(1:9,13,10:12)]
head(malt3)

### add this to the rows with <=1 death
dim(malt[!malt$IntDeaths>1, ])
malt4 <- rbind(malt3, malt[!malt$IntDeaths>1, ])
head(malt4)

table(malt4$IntDeaths, malt4$Censored) # this will check if there are censored flies among the deaths
# eliminate the censored flies in the deaths data:
malt4$Censored <- 0
table(malt4$IntDeaths, malt4$Censored) 
## keep only the rows with a death:
malt4 <- malt4[malt4$IntDeaths == 1, ]
table(malt4$IntDeaths)

### do this same procedure with censored flies######:
### get rows with >1 Censored
malt5 <- malt[malt$Censored>1, ]
head(malt5)
### expand to give one row per censored fly per time
malt6 <- expandRows(malt5, "Censored")
head(malt6)
####### Add back the Censored column, but now all values = 1
malt6$Censored <- 1
head(malt6)
head(malt)
malt6 <- malt6[, c(1:8,13,9:12)]
head(malt6)

### add this to the rows with <=1 Censored
dim(malt4[!malt4$Censored>1, ])
malt7 <- rbind(malt6, malt4[!malt4$Censored>1, ])
head(malt7)

table(malt7$Censored)
# eliminate the deaths from these censored-only data:
malt7$IntDeaths <- 0
table(malt7$Censored, malt7$IntDeaths)
# eliminate rows where censored is not 1
malt7[malt7$Censored == 1, ] -> malt7

# you now have two data frames, one with all of the censored flies, and the other with all of the deaths.  You will rbind these two dfs and this will give a column (IntDeaths) that contains either 1 or 0.  The 1 is when a death occured and 0 is when a censored fly left the experiment.
malt8 <- rbind(malt4, malt7)

#### make a new df$event column ...which should just be the same as IntDeaths because when it is 0 is when a fly was censored.
malt8$event <- malt8$IntDeaths
table(malt8$event)
head(malt8)
malt8 <- separate(malt8, UniqueName, into=c('treatment', 'maltose', 'line'), sep = '/', remove=F)
head(malt8)
str(malt8)
malt8[, 'treatment'] <- as.factor(malt8[, 'treatment'])
malt8[, 'maltose'] <- as.factor(malt8[, 'maltose'])
levels(malt8$maltose) <- c("1.0%", "2.0%", "0.5%", "0%")
malt8$maltose <- factor(malt8$maltose, levels = c("0%", "0.5%", "1.0%", "2.0%"))
head(malt8)
str(malt8)

### get survival perameters (e.g. mean, median, etc) with survfit

malt.stats <- print(survfit(Surv(AgeH, event) ~ malt8$UniqueName, data=malt8), print.rmean=T) # including print((__________), print.rmean=T) gives the mean too  

##### to make the output of survfit into a dataframe:
as.data.frame(summary(malt.stats)$table) -> malt.stats
head(malt.stats)

### make a df of the survfit stats: split $ID into line, sex, treatment, rearrange columns and make some things factors:
head(malt.stats)
malt.stats$ID <- rownames(malt.stats)
malt.stats <- separate(malt.stats, ID, into=c('treatment', 'maltose', 'line'), sep = '/', remove=F)
head(malt.stats)
malt.stats <- separate(malt.stats, treatment, into=c('junk', 'treatment'), sep = '=', remove=T)
head(malt.stats)
malt.stats <- malt.stats[ ,-c(10:11)]
head(malt.stats)
malt.stats <- malt.stats[ ,c(10:12, 1:9)]
head(malt.stats)
rownames(malt.stats)=NULL
colnames(malt.stats) <- c("treatment", "maltose", 'line', "records", "n.max",  "n.start", "event", "rmean", "se.rmean", "median", "0.95LCL", "0.95UCL") ## this step is done to get rid of the astrisks in two of the column names
head(malt.stats)
malt.stats[, 'maltose'] <- as.factor(malt.stats[, 'maltose'])
levels(malt.stats$maltose) <- c("1.0%", "2.0%", "0.5%", "0%")
malt.stats$maltose <- factor(malt.stats$maltose, levels = c("0%", "0.5%", "1.0%", "2.0%"))
malt.stats[, 'treatment'] <- as.factor(malt.stats[, 'treatment'])
str(malt.stats)
table(malt.stats$line)
table(malt.stats$event, malt.stats$treatment)
par(mfrow=c(1,1))
subset(malt.stats, treatment == 'peroxide') -> treats

par(mfrow=c(1,2))
head(treats)
levels(treats$maltose)
plot(rmean ~ maltose, data=treats)
plot(rmean ~ maltose, data=treats[treats$line == 'RAL_161', ]) 
treats$numeric.maltose <- c(1, 1, 1, 1, 2, 2, 2, 2, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0)
head(treats)

## order the maltose levels
order(treats$numeric.maltose)
head(treats)
treats[order(treats[ ,13]) ,] -> treats

par(mfrow=c(1,1))
par(mar=c(6,6,6,6))
plot(rmean ~ numeric.maltose, data=treats[treats$line == 'RAL_161', ], ylim=c(40, 145), col='white', ylab = 'mean lifespan (h)', xlab='supplemental maltose (%)', cex.lab = 2, cex.axis = 1.5)
lines(rmean ~ numeric.maltose, data=treats[treats$line == 'RAL_161',] , type='b', pch=19, lwd=2, col = mycol4[1])
lines(rmean ~ numeric.maltose, data=treats[treats$line == 'RAL_239',] , type='b', pch=19, lwd=2, col = mycol4[2])
lines(rmean ~ numeric.maltose, data=treats[treats$line == 'RAL_320',] , type='b', pch=19, lwd=2, col = mycol4[3])
lines(rmean ~ numeric.maltose, data=treats[treats$line == 'RAL_627',] , type='b', pch=19, lwd=2, col = mycol4[4])

par(xpd=TRUE)
legend(2.2,105 ,c(' RAL 161', ' RAL 239', ' RAL 320', ' RAL 627'), pch = c(19), col=mycol4, bty='n', y.intersp = 1.2)

# Replot, leaving out RAL_161 because we have complete data for the other three in the next experiment (data on RAL_161 in 2% glucose was not obtained - insufficient flies)

subset(treats, line=='RAL_239') -> RAL239
subset(treats, line=='RAL_320') -> RAL320
subset(treats, line=='RAL_627') -> RAL627

par(mfrow=c(1,1))
par(mar=c(6,6,2,6))

plot(rmean ~ numeric.maltose, data=treats[treats$line == 'RAL_161', ], ylim=c(40, 145), col='white', xlab='', ylab='', cex.lab = 2, cex.axis = 1.5, las=1)
arrows(x0 = RAL239$numeric.maltose, y0 = RAL239$rmean, x1 = RAL239$numeric.maltose, y1 = RAL239$rmean + RAL239$se.rmean,  angle = 90, length = 0.05, lwd = 1)
arrows(x0 = RAL239$numeric.maltose, y0 = RAL239$rmean, x1 = RAL239$numeric.maltose, y1 = RAL239$rmean - RAL239$se.rmean,  angle = 90, length = 0.05, lwd = 1)
lines(rmean ~ numeric.maltose, data=treats[treats$line == 'RAL_239',] , type='b', pch=19, lwd=2, col = mycol4[2], cex=1.5)

arrows(x0 = RAL320$numeric.maltose, y0 = RAL320$rmean, x1 = RAL320$numeric.maltose, y1 = RAL320$rmean + RAL320$se.rmean,  angle = 90, length = 0.05, lwd = 1)
arrows(x0 = RAL320$numeric.maltose, y0 = RAL320$rmean, x1 = RAL320$numeric.maltose, y1 = RAL320$rmean - RAL320$se.rmean,  angle = 90, length = 0.05, lwd = 1)
lines(rmean ~ numeric.maltose, data=treats[treats$line == 'RAL_320',] , type='b', pch=19, lwd=2, col = mycol4[3], cex=1.5)

arrows(x0 = RAL627$numeric.maltose, y0 = RAL627$rmean, x1 = RAL627$numeric.maltose, y1 = RAL627$rmean + RAL627$se.rmean,  angle = 90, length = 0.05, lwd = 1)
arrows(x0 = RAL627$numeric.maltose, y0 = RAL627$rmean, x1 = RAL627$numeric.maltose, y1 = RAL627$rmean - RAL627$se.rmean,  angle = 90, length = 0.05, lwd = 1)
lines(rmean ~ numeric.maltose, data=treats[treats$line == 'RAL_627',] , type='b', pch=19, lwd=2, col = mycol4[4], cex=1.5)

par(xpd=TRUE)
legend(2.2,105 ,c(' RAL 239', ' RAL 320', ' RAL 627'), pch = c(19), col=mycol4[2:4], bty='n', y.intersp = 1.2)

set.seed(1)

x <- 1:10
y <- sample(1:15, 10)

par(mar = c(5.1, 2.1, 4.1, 4.1))
plot(x, y, yaxt="n", ylab=NA)
axis(4)


############

### test for the effect of different levels of maltose:
aov(rmean ~ maltose + line, data=treats) -> aov.malt.line
summary(aov.malt.line)
TukeyHSD(x=aov.malt.line, 'maltose')  #get pairwise tests between the maltose levels along with adjusted pValues
plot(aov.malt.line) # hit enter in the console below to see some diagnostic plots for the linear model

################################################################################################




