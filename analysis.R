##No probe first
no_probe <- read.csv("noprobe.csv", header = T)

#calculate the mRNA per Area
no_probe$mRNA <- no_probe$Counts / no_probe$Area

#take the average of the regions for each spot of each patient and generate a
#data frame that can be used for the rest of the analysis

no_probe_average <<- data.frame() #output data frame

for (i in unique(no_probe$Sample.ID)){
  Sample.ID <<- i
  for (j in unique(no_probe$Spot.ID)){
    Spot.ID <<- j
    calcdf <- subset(no_probe, no_probe$Sample.ID == i & no_probe$Spot.ID == j)
    if (nrow(calcdf) > 1){ # removes a sample if there is less than 2 regions
      Average <- mean(calcdf$mRNA)
      Stage <- unique(calcdf$Stage)
      new <- data.frame(Sample.ID, Stage, Spot.ID, Average)
      no_probe_average <<- rbind(no_probe_average, new) #generating the new data frame
    }
  }
}

##Now for the probe

probe <- read.csv("pgc1beta.csv", header = T)

#calculate the mRNA per Area
probe$mRNA <- probe$Counts / probe$Area

#take the average of the regions for each spot of each patient and generate a
#data frame that can be used for the rest of the analysis

probe_average <<- data.frame() #output data frame

for (i in unique(probe$Sample.ID)){
  Sample.ID <<- i
  for (j in unique(probe$Spot.ID)){
    Spot.ID <<- j
    calcdf <- subset(probe, probe$Sample.ID == i & probe$Spot.ID == j)
    if (nrow(calcdf) > 1){ # removes a sample if there is less than 2 regions
      Average <- mean(calcdf$mRNA)
      Stage <- unique(calcdf$Stage)
      new <- data.frame(Sample.ID, Stage, Spot.ID, Average)
      probe_average <<- rbind(probe_average, new) #generating the new data frame
    }
  }
}

##Background Correction
#Match patient and subtract no probe average from probe if available

probe_nobkgd <<- data.frame() #output data frame

for (i in unique(probe_average$Sample.ID)){
  Sample.ID <<- i
  for (j in unique(probe_average$Spot.ID)){
    Spot.ID <<- j
    calcdfp <- subset(probe_average, probe_average$Sample.ID == i & probe_average$Spot.ID == j)
    calcdfnp <- subset(no_probe_average, no_probe_average$Sample.ID == i & no_probe_average$Spot.ID == j)
    if (nrow(calcdfnp) > 0 & nrow(calcdfp) > 0){ #checks to make sure there are values. if not, no match
        pmnp <- calcdfp$Average - calcdfnp$Average #probe minus no probe
        Stage <- unique(calcdfp$Stage)
        new <- data.frame(Sample.ID, Stage, Spot.ID, pmnp)
        probe_nobkgd <<- rbind(probe_nobkgd, new)
    }
  }
}

##Compute Fold Change

stage_average <- tapply(probe_nobkgd$pmnp, list(probe_nobkgd$Stage), mean, na.rm = T)

#create new data frame containing old info plus fold change
probe_fc <- probe_nobkgd
probe_fc$fc <- probe_nobkgd$pmnp / stage_average["Normal"]

stage_average_fc <- tapply(probe_fc$fc, list(probe_fc$Stage), mean, na.rm = T)

##output fold change data
write.csv(probe_fc, file = "pgc1beta-fc.csv")

##stats - anova
#test for normality
shapiro.test(probe_fc$fc)
#One-way ANOVA and Tukey post-hoc comparison test
probe_aov <- aov(fc~Stage, probe_fc)
summary(probe_aov)
TukeyHSD(probe_aov)

##graph time

plot(probe_fc$Stage, probe_fc$fc)

