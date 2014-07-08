#PGC1beta plots

#packages needed to make pretty plots
library(ggplot2, plyr)

#raw data
ourData <- subset(probe, probe$TMA != "TMA3")
wangData <- subset(probe, probe$TMA == "TMA3")

#background subtracted and averaged to sample
ourDataNBG <- subset(probe_PA, probe_PA$TMA != "TMA3")
wangDataNBG <- subset(probe_PA, probe_PA$TMA == "TMA3")

#some threshold info
threshold_TMA <- tapply(probe$Threshold, list(probe$TMA), mean, na.rm = T)
threshold_ourData <-tapply(ourData$Threshold, list(ourData$Stage), mean, na.rm = T)
threshold_wangData <-tapply(wangData$Threshold, list(wangData$Stage), mean, na.rm = T)

#Stats on our data TMA's 1,2, 4,5,6  NO 3
ourDataNBG_aov <- aov(Average~Stage, ourDataNBG)
summary(ourDataNBG_aov)
TukeyHSD(ourDataNBG_aov)

#Subsetting our data wihtout TMA into separte stages to run t tests on each
ourDataNBG_Normal <- subset(ourDataNBG, ourDataNBG$Stage == "Normal")
ourDataNBG_Stage1 <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 1")
ourDataNBG_Stage2 <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 2")
ourDataNBG_Stage3 <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 3")
ourDataNBG_Stage4 <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 4")
ourDataNBG_Stage3Met <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 3 Met")
ourDataNBG_Stage4Met <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 4 Met")
ourDataNBG_TA <- subset(ourDataNBG, ourDataNBG$Stage == "Tubular Adenoma")


#the t tests
ourDataNBG_ttest <- t.test(ourDataNBG_Normal$Average, ourDataNBG_Stage1$Average)
ourDataNBG_ttest2 <- t.test(ourDataNBG_Normal$Average, ourDataNBG_Stage2$Average)
ourDataNBG_ttest3 <- t.test(ourDataNBG_Normal$Average, ourDataNBG_Stage3$Average)
ourDataNBG_ttest4 <- t.test(ourDataNBG_Normal$Average, ourDataNBG_Stage4$Average)
ourDataNBG_ttest5 <- t.test(ourDataNBG_Normal$Average, ourDataNBG_Stage3Met$Average)
ourDataNBG_ttest6 <- t.test(ourDataNBG_Normal$Average, ourDataNBG_Stage4Met$Average)
ourDataNBG_ttest7 <- t.test(ourDataNBG_Normal$Average, ourDataNBG_TA$Average)


#data for error bars for plot 1
ourDataNBG_stat <- ddply(ourDataNBG, .(Stage), summarise, mean = mean(Average), sd = sd(Average))

#plot 1 - scatter plot with error bars, all data except TMA 3
ggplot(ourDataNBG, aes(x = Stage, y = Average)) +
  geom_point(aes(color = factor(ourDataNBG$TMA)), size = 5) + 
  geom_point(data = ourDataNBG_stat, aes(x = Stage, y = mean), colour = 'black', size = 5) +
  geom_errorbar(data = ourDataNBG_stat, aes(x = Stage, y = Average, ymin = mean - sd, ymax = mean + sd), width = 0.4) +
  labs(title = "PGC1beta in Human Colon Cancer Samples", x = "", y = "mRNA / Area(pixels)", colour = "")

ggsave(
  "PGC1beta plot 1.png",
  dpi = 500
)

#subsetting normals, stage 4, stage 4 mets, no tma 3 samples
ourDataNPM <- subset(ourDataNBG, ourDataNBG$Stage == "Normal" | ourDataNBG$Stage == "Stage 4" | ourDataNBG$Stage == "Stage 4 Met")

#stats for plot 2
ourDataNPM_aov <- aov(Average~Stage, ourDataNPM)
summary(ourDataNPM_aov)
TukeyHSD(ourDataNPM_aov)

#for error bars for plot 2
ourDataNPM_stat <- ddply(ourDataNPM, .(Stage), summarise, mean = mean(Average), sd = sd(Average))

#plot 2 - Normals vs Stage 4 Primrary vs Stage 4 Met, no TMA 3
ggplot(ourDataNPM, aes(x = Stage, y = Average)) +
  geom_point(size = 5) + 
  geom_point(data = ourDataNPM_stat, aes(x = Stage, y = mean), colour = 'red', size = 6) +
  geom_errorbar(data = ourDataNPM_stat, aes(x = Stage, y = Average, ymin = mean - sd, ymax = mean + sd), width = 0.4) +
  labs(title = "PGC1beta in Human Colon Cancer Samples", x = "", y = "mRNA / Area(pixels)", colour = "")

ggsave(
  "PGC1beta plot 2.png",
  width = 5.0,
  height = 5.0,
  dpi = 500
)

#error bar info for plot 3
probe_PA_stat <- ddply(probe_PA, .(Stage), summarise, mean = mean(Average), sd = sd(Average))

#plot of all data, including TMA 3
ggplot(probe_PA, aes(x = Stage, y = Average)) +
  geom_point(aes(color = factor(TMA)), size = 5, position = "dodge") + 
  geom_point(data = probe_PA_stat, aes(x = Stage, y = mean), colour = 'black', size = 6) +
  geom_errorbar(data = probe_PA_stat, aes(x = Stage, y = Average, ymin = mean - sd, ymax = mean + sd), width = 0.4) +
  labs(title = "PGC1beta in Human Colon Cancer Samples", x = "", y = "mRNA / Area(pixels)", colour = "")

ggsave(
  "PGC1beta plot 3.png",
  dpi = 500
)