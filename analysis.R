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
    if (nrow(calcdfnp) > 0 & nrow(calcdfp) > 0){
        pmnp <- calcdfp$Average - calcdfnp$Average #probe minus no probe
        Stage <- unique(calcdfp$Stage)
        new <- data.frame(Sample.ID, Stage, Spot.ID, pmnp)
        probe_nobkgd <<- rbind(probe_nobkgd, new)
    }
  }
}

##Take average of each patient

##Compute fold change based on normal