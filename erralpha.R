#ERRalpha plots

jpeg("ERRalphaAllFoldChange.jpg")
plot(probe_fc$Stage, probe_fc$fc, main = "ERRalpha compared to Normal")
dev.off()

jpeg("ERRalphaAllFoldChangeBar.jpg")
barplot(stage_average_fc)
dev.off()

jpeg("RawData.jpg")
plot(probe$Stage, probe$mRNA, main = "Raw Data, no background subtraction")
dev.off()

threshold_probe <- tapply(probe$Threshold, list(probe$Stage), mean, na.rm = T)
ourData <- subset(probe, probe$TMA != "TMA3")
wangData <- subset(probe, probe$TMA == "TMA3")

jpeg("NoTMA3Raw.jpg")
plot(ourData$Stage, ourData$mRNA, main = "Raw Data, no TMA3")
dev.off()

jpeg("TMA3Raw.jpg")
plot(wangData$Stage, wangData$mRNA, main = "Wang TMA only")
dev.off()

jpeg("ByTMA.jpg")
plot(probe_fc$TMA, probe_fc$fc, main = "By TMA")
dev.off()

jpeg("BackgroundSubAll.jpg")
plot(probe_nobkgd$Stage, probe_nobkgd$pmnp, main = "Background Subtracted All")
dev.off()

ourDataNBG <- subset(probe_nobkgd, probe_nobkgd$TMA != "TMA3")
wangDataNBG <- subset(probe_nobkgd, probe_nobkgd$TMA == "TMA3")

jpeg("NoTMA3NBKG.jpg")
plot(ourDataNBG$Stage, ourDataNBG$pmnp, main = "Background Subtracted, no TMA3")
dev.off()

jpeg("TMA3NBKG.jpg")
plot(wangDataNBG$Stage, wangDataNBG$pmnp, main = "Background Stubtacted, Wang TMA only")
dev.off()

threshold_TMA <- tapply(probe$Threshold, list(probe$TMA), mean, na.rm = T)
threshold_ourData <-tapply(ourData$Threshold, list(ourData$Stage), mean, na.rm = T)
threshold_wangData <-tapply(wangData$Threshold, list(wangData$Stage), mean, na.rm = T)

#Stats on our data TMA's 1,2, 4,5,6  NO 3
ourDataNBG_aov <- aov(pmnp~Stage, ourDataNBG)
summary(ourDataNBG_aov)
TukeyHSD(ourDataNBG_aov)

ourDataNBG_Normal <- subset(ourDataNBG, ourDataNBG$Stage == "Normal")
ourDataNBG_Stage1 <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 1")
ourDataNBG_Stage2 <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 2")
ourDataNBG_Stage3 <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 3")
ourDataNBG_Stage4 <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 4")
ourDataNBG_Stage3Met <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 3 Met")
ourDataNBG_Stage4Met <- subset(ourDataNBG, ourDataNBG$Stage == "Stage 4 Met")
ourDataNBG_TA <- subset(ourDataNBG, ourDataNBG$Stage == "Tubular Adenoma")



ourDataNBG_ttest <- t.test(ourDataNBG_Normal$pmnp, ourDataNBG_Stage1$pmnp)
ourDataNBG_ttest2 <- t.test(ourDataNBG_Normal$pmnp, ourDataNBG_Stage2$pmnp)
ourDataNBG_ttest3 <- t.test(ourDataNBG_Normal$pmnp, ourDataNBG_Stage3$pmnp)
ourDataNBG_ttest4 <- t.test(ourDataNBG_Normal$pmnp, ourDataNBG_Stage4$pmnp)
ourDataNBG_ttest5 <- t.test(ourDataNBG_Normal$pmnp, ourDataNBG_Stage3Met$pmnp)
ourDataNBG_ttest6 <- t.test(ourDataNBG_Normal$pmnp, ourDataNBG_Stage4Met$pmnp)
ourDataNBG_ttest7 <- t.test(ourDataNBG_Normal$pmnp, ourDataNBG_TA$pmnp)


library(ggplot2)

ggplot(ourDataNBG, aes(x = ourDataNBG$Stage, y = ourDataNBG$pmnp)) +
  geom_point(aes(color = factor(ourDataNBG$TMA))) + 
  labs(title = "ERRalpha in Human Colon Cancer Samples", x = "", y = "mRNA / Area(pixels)", colour = "")

ourDataNPM <- subset(ourDataNBG, ourDataNBG$Stage == "Normal" | ourDataNBG$Stage == "Stage 4" | ourDataNBG$Stage == "Stage 4 Met")

ourDataNPM_stat <- ddply(ourDataNPM, .(Stage), summarise, mean = mean(pmnp), sd = sd(pmnp))

ggplot(ourDataNPM, aes(x = Stage, y = pmnp)) +
  geom_point(size = 3) + 
  geom_point(data = ourDataNPM_stat, aes(x = Stage, y = mean), colour = 'red', size = 5) +
  geom_errorbar(data = ourDataNPM_stat, aes(x = Stage, y = pmnp, ymin = mean - sd, ymax = mean + sd), width = 0.4) +
  labs(title = "ERRalpha in Human Colon Cancer Samples", x = "", y = "mRNA / Area(pixels)", colour = "")

ggsave(
  "ERRalpha plot 2.png",
  width = 5.0,
  height = 5.0,
  dpi = 500
)

probe_nobkgd_stat <- ddply(probe_nobkgd, .(Stage), summarise, mean = mean(pmnp), sd = sd(pmnp))

ggplot(probe_nobkgd, aes(x = Stage, y = pmnp)) +
  geom_point(aes(color = factor(TMA)), size = 5, position = "dodge") + 
  geom_point(data = probe_nobkgd_stat, aes(x = Stage, y = mean), colour = 'black', size = 6) +
  geom_errorbar(data = probe_nobkgd_stat, aes(x = Stage, y = pmnp, ymin = mean - sd, ymax = mean + sd), width = 0.4) +
  labs(title = "ERRalpha in Human Colon Cancer Samples", x = "", y = "mRNA / Area(pixels)", colour = "")

ggsave(
  "ERRalpha plot 3.png",
  dpi = 500
)

ggplot(no_probe, aes(x = TMA, y = mRNA)) +
  geom_boxplot() + 
  labs(title = "No Probe raw values by TMA")

ggsave(
  "No Probe Raw Values.png",
  dpi = 500
)
