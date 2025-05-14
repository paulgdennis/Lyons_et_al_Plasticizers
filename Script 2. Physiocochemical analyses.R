### SCRIPT 2. Lyons et al., (2025) 

# Paul Dennis (p.dennis@uq.edu.au) and Rebecca Lyons (r.lyons@uq.edu.au)

### List of scripts ###
# Script 1: Data input and sub-setting
# Script 2: Physicochemical analyses
# Script 3: Biomass, alpha, beta diversity
# Script 4: Indicator analyses
# Script 5: Differential abundance testing

############################
# 2 Physicochemical analyses

source('Script 1. Data input and subsetting.R')
library(sciplot)

### GPC (Mn, Mw, PD)

# Plot data from 0 and 14 days as bars with SEM and individual data points

# Identify where it where the data points are on the x axis. 
nTimes <- 2
chem_gpc$i_Plasticiser <- as.numeric(factor(chem_gpc$Plasticiser , 
                                            levels = c('Ctrl','TBAC','TEC','DBP')))
chem_gpc$i_Time <- as.numeric(factor(chem_gpc$Time))
chem_gpc$x_barplot <- (chem_gpc$i_Plasticiser - 1) * (nTimes + 0.1 + 0.2) + (chem_gpc$i_Time - 1) + 0.5 + 0.2 + 0.1
chem_gpc$x_barplot[chem_gpc$i_Time == 1] <- chem_gpc$x_barplot[chem_gpc$i_Time == 1] - 0.1

# Bargraph with data points overlaid (substitute variable for Mw, dispersity)
bargraph.CI(x.factor = Plasticiser, 
            response = Mn, 
            group = Time,
            data = chem_gpc, legend = FALSE, err.width = 0.05, 
            col = c('#339e2bff','#339e2bff','#b0de8aff','#b0de8aff','#a6cce3ff','#a6cce3ff','#377EB8','#377EB8'), 
            space = c(0.1, 0.2) , 
            ylim = c(0 , 1.05 * max(chem_gpc$Mn))) 
points(chem_gpc$x_barplot , chem_gpc$Mn , pch = 16 , cex = 0.4 , col = 'black')

### NMR

# Plot NMR data from 0 and 14 days as bars with SEM and individual data points

##Identify where it where the datapoints are on the x axis. 
nTimes <- 2
chem_nmr$i_Plasticiser <- as.numeric(factor(chem_nmr$Plasticiser, 
                                            levels = c('Ctrl','TBAC','TEC','DBP')))
chem_nmr$i_Time <- as.numeric(factor(chem_nmr$Time))
chem_nmr$x_barplot <- (chem_nmr$i_Plasticiser - 1) * (nTimes + 0.1 + 0.2) + (chem_nmr$i_Time - 1) + 0.5 + 0.2 + 0.1
chem_nmr$x_barplot[chem_nmr$i_Time == 1] <- chem_nmr$x_barplot[chem_nmr$i_Time == 1] - 0.1

#Make bargraph and overlay datapoints
bargraph.CI(chem_nmr$Plasticiser, chem_nmr$plasticiser_pc,chem_nmr$Time,
            legend = FALSE, xlab = "", ylab = "Plasticiser content", err.width = 0.05, 
            col=c('#339e2bff','#339e2bff','#b0de8aff','#b0de8aff','#a6cce3ff','#a6cce3ff','#377EB8','#377EB8'), 
            space = c(0.1, 0.2), ylim = c(0,1.02 * max(chem_nmr$plasticiser_pc)))
points(chem_nmr$x_barplot , chem_nmr$plasticiser_pc , pch = 16 , cex = 0.4 , col = 'black')

### Two sample t-test between 0 and 14 days

# Subset into before and after
# Substitute Mn with Mw or dispersity, substitue dataframe with chem_nmr and plasticiser_pc for nmr data.
chem_gpc.before <- chem_gpc[chem_gpc$Time == "before",]
chem_gpc.after <- chem_gpc[chem_gpc$Time == "after",]

group1 <- chem_gpc.before$Mn[which(chem_gpc.before$Plasticiser == "Ctrl")]
group2 <- chem_gpc.after$Mn[which(chem_gpc.after$Plasticiser == "Ctrl")]
t.test(group1, group2)

group1 <- chem_gpc.before$Mn[which(chem_gpc.before$Plasticiser == "TBAC")]
group2 <- chem_gpc.after$Mn[which(chem_gpc.after$Plasticiser == "TBAC")]
t.test(group1, group1)

group1 <- chem_gpc.before$Mn[which(chem_gpc.before$Plasticiser == "TEC")]
group2 <- chem_gpc.after$Mn[which(chem_gpc.after$Plasticiser == "TEC")]
t.test(group1, group2)

group1 <- chem_gpc.before$Mn[which(chem_gpc.before$Plasticiser == "DBP")]
group2 <- chem_gpc.after$Mn[which(chem_gpc.after$Plasticiser == "DBP")]
t.test(group1, group2)
