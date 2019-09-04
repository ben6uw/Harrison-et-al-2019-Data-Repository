mycol <- c("#8B8378", "#EE9A00") 
mycol2 <- c(mycol,"#483D8B")
mycol4 = c("#EEC900", "#CD6600", "#87CEFA", "#104E8B")
mycol3 <- mycol4

setwd("/Users/ben/Google Drive/Documents/Peroxide/Elises.Feeding.Experiments")
list.files()

perfeed2 <- read.csv('peroxide.feeding.assay.csv', header=T, stringsAsFactors=T)
names(perfeed2)
perfeed2$linefood <- paste(perfeed2$line, perfeed2$food, sep=".")
names(perfeed2) <- c("vial", "line", "food",  "height.consumed", "alive", "linefood")
perfeed2$volume.consumed = perfeed2$height.consumed * 0.441786467
head(perfeed2)
str(perfeed2)
ifelse(is.na(perfeed2$line), 'no flies', perfeed2$line) -> perfeed2$line
as.factor(paste(perfeed2$line, perfeed2$food)) -> perfeed2$linefood
par(mar=c(12,6,6,6))
par(mfrow=c(1,1))

plot(perfeed2$volume.consumed ~ as.factor(perfeed2$linefood), main= "Peroxide Feeding Assay", ylab='Food Consumed (μL)', xlab=NULL, cex.main=2, cex.lab= 2, cex.axis=2, outline=F, las=2, col= c(mycol4[c(2,2,3,3,4,4)],0,0), boxlwd = 1.5)

box(lwd=1.5)

stripchart(perfeed2$volume.consumed ~ perfeed2$linefood,vertical = T, add= T, method = 'jitter', cex=1, pch=c(19, 21), bg='white')


levels(perfeed2$linefood)

wilcox.test(perfeed2$volume.consumed[perfeed2$linefood == '239 peroxide'], perfeed2$volume.consumed[perfeed2$linefood == '239 water'])
wilcox.test(perfeed2$volume.consumed[perfeed2$linefood == '320 peroxide'], perfeed2$volume.consumed[perfeed2$linefood == '320 water'])
wilcox.test(perfeed2$volume.consumed[perfeed2$linefood == '627 peroxide'], perfeed2$volume.consumed[perfeed2$linefood == '627 water'])
wilcox.test(perfeed2$volume.consumed[perfeed2$linefood == '239 peroxide'], perfeed2$volume.consumed[perfeed2$linefood == 'no flies peroxide'])

head(perfeed2)
table(perfeed2$linefood)

feed3 <- aggregate(perfeed2$volume.consumed,
                   by = list(perfeed2$linefood),
                   FUN = function(x) c(mean = mean(x), sd = sd(x)))

data.frame(feed3$x) -> feed3
feed3
feed3[c(2,1,4,3,6,5,8,7), ] -> feed4
feed4

max(feed4$mean) + max(feed4$sd) -> tippytop

par(mar=c(6,8,6,6))
plot <- barplot(height = feed4$mean[1:7], col=c(mycol4[4], 0, mycol4[4], 0, mycol4[4], 0, 'grey'), lwd=1, las=1, space=c(0.4,0,0.4,0), border=T, ylim = c(0, 1.1*tippytop), ann=F, axes=F)
axis(2,at=seq(0,1.1*tippytop,1), line=1.2, cex.axis=1.5, las=2)
arrows(plot, feed4$mean[1:7], y1 = feed4$mean[1:7] + feed4$sd[1:7], lwd = 1.5, angle=90, code=2, length=0.05)

title(ylab = "food consumed (uL)", line=4.5, cex.lab=2)
title(xlab = "    Ral_239  Ral_320  Ral_627 no flies", cex.lab=2, line=2)

