##No probe first
no_probe <- read.csv("noprobe.csv", header = T)

#pixels to um - 32 px = 0.01mm  
no_probe$AreaUM <- no_probe$Area / (32 * 32 / 10)

#calculate the mRNA per Area
no_probe$mRNA <- no_probe$Count / no_probe$AreaUM

#take the average of the regions for each spot of each patient and generate a
#data frame that can be used for the rest of the analysis

no_probe_average <<- data.frame() #output data frame

for (k in unique(no_probe$TMA)){
  TMA <<- k
  for (i in unique(no_probe$Sample)){
    Sample <<- i
    for (j in unique(no_probe$Spot)){
      Spot <<- j
      calcdf <- subset(no_probe, no_probe$TMA == k & no_probe$Sample == i & no_probe$Spot == j)
      if (nrow(calcdf) > 1){ # removes the sample if there is less than 2 regions
        Average <- mean(calcdf$mRNA)
        Stage <- unique(calcdf$Stage)
        TMA <- unique(calcdf$TMA)
        tmp <- data.frame(Sample, TMA, Stage, Spot, Average)
        no_probe_average <<- rbind(no_probe_average, tmp) #generating the new data frame
      }
    }
  }
}

##Now for the probe

probe <- read.csv("pgc1beta.csv", header = T)

#pixels to um - 32 px = 0.01mm  
probe$AreaUM <- probe$Area / (32 * 32 / 10)

#calculate the mRNA per Area
probe$mRNA <- probe$Count / probe$AreaUM

#take the average of the regions for each spot of each patient and generate a
#data frame that can be used for the rest of the analysis

probe_average <<- data.frame() #output data frame

for (k in unique(probe$TMA)){
  TMA <<- k
  for (i in unique(probe$Sample)){
    Sample <<- i
    for (j in unique(probe$Spot)){
      Spot <<- j
      calcdf <- subset(probe, probe$TMA == k & probe$Sample == i & probe$Spot == j)
      if (nrow(calcdf) > 1){ # removes a sample if there is less than 2 regions
        Average <- mean(calcdf$mRNA)
        Stage <- unique(calcdf$Stage)
        TMA <- unique(calcdf$TMA)
        tmp <- data.frame(Sample, TMA, Stage, Spot, Average)
        probe_average <<- rbind(probe_average, tmp) #generating the new data frame
      }
    }
  }
}

##Background Correction
#Match patient and subtract no probe average from probe if available

probe_nobkgd <<- data.frame() #output data frame

for (k in unique(probe_average$TMA)){
  TMA <<- k
  for (i in unique(probe_average$Sample)){
    Sample <<- i
    for (j in unique(probe_average$Spot)){
      Spot <<- j
      calcdfp <- subset(probe_average, probe_average$TMA == k & probe_average$Sample == i & probe_average$Spot == j)
      calcdfnp <- subset(no_probe_average, no_probe_average$TMA == k & no_probe_average$Sample == i & no_probe_average$Spot == j)
      if (nrow(calcdfnp) > 0 & nrow(calcdfp) > 0){ #checks to make sure there are values. if not, no match
        pmnp <- calcdfp$Average - calcdfnp$Average #probe minus no probe
        Stage <- unique(calcdfp$Stage)
        TMA <- unique(calcdfp$TMA)
        tmp <- data.frame(Sample, TMA, Stage, Spot, pmnp)
        probe_nobkgd <<- rbind(probe_nobkgd, tmp)
      }
    }
  }
}

##Average each patient
probe_PA <<- data.frame() #Patient Average

for (k in unique(probe_nobkgd$TMA)){
  TMA <<- k
  for (i in unique(probe_nobkgd$Sample)){
    Sample <<- i
    calcdf <- subset(probe_nobkgd, probe_nobkgd$TMA == k & probe_nobkgd$Sample == i)
    if (nrow(calcdf) > 0) {
      Average <- mean(calcdf$pmnp)
      Stage <- unique(calcdf$Stage)
      tmp <- data.frame(TMA, Sample, Stage, Average)
      probe_PA <<- rbind(probe_PA, tmp) #generating the new data frame
    }
  }
}




##Compute Fold Change

stage_average <- tapply(probe_PA$Average, list(probe_PA$Stage), mean, na.rm = T)

#create new data frame containing old info plus fold change
probe_fc <- probe_PA
probe_fc$fc <- probe_PA$Average / stage_average["Normal"]

stage_average_fc <- tapply(probe_fc$fc, list(probe_fc$Stage), mean, na.rm = T)

##output fold change data
#write.csv(probe_fc, file = "pgc1beta-fc.csv")

##stats - anova
#test for normality
shapiro.test(probe_fc$fc)
#One-way ANOVA and Tukey post-hoc comparison test
probe_aov <- aov(fc~Stage, probe_fc)
summary(probe_aov)
TukeyHSD(probe_aov)

