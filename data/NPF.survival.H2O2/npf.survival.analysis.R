# analysis of the NPF mutant experiment(s) that Brian Y. Chung conducted in Scott Pletcher's lab Spring 2019:
library(splitstackshape)
library(survival)

read.table("../SurvData_2019-04-12 0820.csv", header=T, sep=',', stringsAsFactors = F) -> dat 

head(dat)
unique(dat$strata)
cSplit(dat, 'strata', sep='=', drop=F) -> dat
dat[ ,c(1, 10, 2:8)] -> dat
head(dat)
table(dat$n.event) # rows are deaths, most are single deaths, but four rows record two deaths
expandRows(dat, 'n.event')-> dat
dat$event <- 1 # makes an 'event' column, where each 'event' is a death

dat$sex <- NA
dat$sex[grep('F', dat$strata_2)] <- 'female'
dat$sex[grep('M', dat$strata_2)] <- 'male'
dat$treatment <- 'water'
dat$treatment[grep('H', dat$strata_2)] <- 'peroxide'
dat$genotype <- 'CantonS'
dat$genotype[grep('SK', dat$strata_2)] <- 'NPF-SK1'
table(dat$genotype)
table(dat$treatment, dat$genotype)
head(dat)
dat[dat$treatment == 'peroxide', ] -> npf

head(npf)
table(npf$event, npf$strata_2)
table(npf$n.event)
hist(log(npf$time), 10, col=8) # histogram of log(time of deaths)

## some fun wth plots:
library(lattice)
bwplot(time ~ genotype,  data = subset(npf, sex == 'female'), layout = c(3, 1), panel = panel.violin)
densityplot(~ time, groups = genotype, data = subset(npf, sex == 'female'), plot.points = FALSE, auto.key = TRUE)
#########

plot(survfit(Surv(time, event) ~ genotype, data = subset(npf, sex == 'female'), conf.type='none'), col=0, cex.main=2, las=1)
lines(survfit(Surv(time, event) ~ genotype, data = subset(npf, sex == 'female'), conf.type='none'), col=1:2, type='l', lwd=2)
