##No probe first

no_probe <- read.csv("noprobe.csv", header = T)

#calculate the mRNA per Area
no_probe$mRNA <- no_probe$Counts / no_probe$Area

#take the average of the regions for each spot of each patient
