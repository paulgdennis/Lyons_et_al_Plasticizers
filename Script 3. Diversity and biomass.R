### SCRIPT 3. Lyons et al., (2025) 

# Paul Dennis (p.dennis@uq.edu.au) and Rebecca Lyons (r.lyons@uq.edu.au)

### List of scripts ###
# Script 1: Data input and sub-setting
# Script 2: Physicochemical analyses
# Script 3: Biomass, alpha, beta diversity
# Script 4: Indicator analyses
# Script 5: Differential abundance testing

# Biomass, alpha, beta diversity, Figures 4, 5, 6, 7, S2, Tables 1, S2, S3)

source('Script 1. Data input and subsetting.R')
library(lsmeans)
library(multcomp)
library(sciplot)
library(vegan)

#######
# 16S #
#######

### Biomass - qPCR

# ANOVA with posthoc 
summary(aov(log2(qPCR$X16s.copies.g.material.1in100) ~ Plasticiser * Compartment, data = qPCR))
cld(lsmeans(
  aov(log2(qPCR$X16s.copies.g.material.1in100) ~ Treat.order * Comp.order, data = env_16S), 
            ~ Treat.order * Comp.order), 
            Letters = letters, 
            adjust = "Tukey")

# Plot data from different treatment combinations as bars with 95% CI and individual data points
nTreatments <- 4
qPCR$i_Comp <- as.numeric(factor(qPCR$Compartment , levels = c('BS','AS','P')))
qPCR$i_Plasticiser <- as.numeric(factor(qPCR$Plasticiser , levels = c('Ctrl' , 'TBAC' , 'TEC' , 'DBP')))
qPCR$x_barplot <- (qPCR$i_Comp - 1) * (nTreatments + 1) + qPCR$i_Plasticiser + 0.5

bargraph.CI(env_16S$Comp.order, log2(qPCR$X16s.copies.g.material.1in100), env_16S$Treat.order,
            legend = FALSE, ylim = c(22,36), xlab = "", ylab = "16S copy number per g material",
            main = "16S copy number per g material", err.width = 0.05, col=plas.cols)
points(qPCR$x_barplot , log2(qPCR$X16s.copies.g.material.1in100), pch = 16 , cex = 0.4 , col = 'black')

### Alpha diversity

# ANOVA with Tukey HSD 
for(i in colnames(env_16S[,5:8])){
  print(i)
  print(summary(aov(env_16S[,i] ~ Treatment * Compartment, data = env_16S)))
}

for(i in colnames(env_16S[,5:8])) {
  print(i)
  print(cld(lsmeans(
    aov(env_16S[,i] ~ Treatment * Compartment, data = env_16S), 
    ~ Compartment * Treatment), 
    Letters = letters), 
    adjust = "Tukey")
}

nTreatments <- 4
env_16S$i_Comp.order <- as.numeric(factor(env_16S$Comp.order, 
                                          levels = c('a.BS','b.AS','c.P')))
env_16S$i_Treatment <- as.numeric(factor(env_16S$Treatment, 
                                         levels = c('Ctrl','TBAC','TEC','DBP')))
env_16S$x_barplot <- (env_16S$i_Comp.order - 1) * (nTreatments + 1) + env_16S$i_Treatment + 0.5
  
for(i in colnames(env_16S[,5:8])){
  bargraph.CI(env_16S$Comp.order, env_16S[,i],env_16S$Treat.order,
              legend = FALSE, xlab = "", ylab = i, main = i, err.width = 0.05, 
              col = plas.cols , ylim = c(0 , 1.4 * max(env_16S[,i])))
  points(env_16S$x_barplot , env_16S[,i] , pch = 16 , cex = 0.4 , col = 'black')
}

### Beta diversity

# PERMANOVA of Hellinger transformed OTUs
adonis2(sqrt(otu_16S) ~ Treatment * Compartment, data = env_16S, method = 'euc', permu = 4999)

# PERMOANOVA of Weighted UniFrac
adonis2(w_unifrac_16S ~ Treatment * Compartment, data = env_16S, method = 'euc', permu = 4999)

### Visualise compositional differences with RDA

# Generate an RDA object constrained by Plasticiser:Compartment
otu_16S.rda <- rda(sqrt(otu_16S) ~ TreatComp, data = env_16S)
ord = otu_16S.rda

# Set colors, size and transparency for site symbols in RDA
addTrans <- function(col, alpha = 1) {
  rgb.mat <- col2rgb(col) / 255
  apply(rgb.mat, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}

compartment.cex <- c(2, 4.5, 8)
cex = compartment.cex[factor(env_16S$Comp.order)]

axis.percent <- function(ordination){
  round((100*eigenvals(ordination)[1:2]/ordination$tot.chi[[1]]),digits=2)
}

# Plot the RDA
scaling.val = 3
plot(ord, 
     type='n', scaling=scaling.val, 
     xlab=paste("RDA1 (",axis.percent(
       ord)[[1]],"%)",sep=""),
     ylab=paste("RDA2 (",axis.percent(
       ord)[[2]],"%)",sep=""))

points(ord, dis='sp', pch=4, col='grey', cex=0.6, scaling=scaling.val)

points(ord, dis='sites', pch=21, scaling = scaling.val,
       bg = addTrans(plas.cols[factor(env_16S$Treat.order)], 180/255),
       cex = compartment.cex[factor(env_16S$Comp.order)])

sd.val = 6
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1])),],
       pch=4,col='darkred',cex=0.6)	
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1]))),],
       pch=4,col='darkred',cex=0.6) 
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2])),],
       pch=4,col='darkred',cex=0.6)
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2]))),],
       pch=4,col='darkred',cex=0.6) 
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1])),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1]))),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2])),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2]))),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)

legend("topleft",legend=unique(factor(env_16S$Treatment)),pch=19,
       col=plas.cols[unique(factor(env_16S$Treat.order))])

# Alternative plot for supplementary figures with different shapes for each compartment:
points(ord, dis='sites', 
       pch = c(21, 22, 24)[factor(env_16S$Comp.order)],  # Different shapes
       scaling = scaling.val,
       bg = addTrans(plas.cols[factor(env_16S$Treat.order)], 180/255),
       cex = 1.5)  # Keep size consistent

##############################
# Generate a heatmap for 16S #
##############################
            
# Calculate the mean relative abundance of each OTUs in each treatment combination
otu_16S_means.tmp <- aggregate(otu_16S, by = list(env_16S$TreatComp), mean)
  row.names(otu_16S_means.tmp) <- otu_16S_means.tmp[,1]
  otu_16S_means.tmp <- otu_16S_means.tmp[,-1]

# Subset the OTU table to populations present at >=X% relative abundance
otus_for_heatmap <- colnames(otu_16S_means.tmp[,which(apply(otu_16S_means.tmp,2,max) >= 0.01)])
hm.tmp <- otu_16S[, otus_for_heatmap] 

# Get the taxonomy
hm.otu.tax <- taxonomy_16S[colnames(hm.tmp),]

hm.otu.tax$full.id <- paste("[",
                            hm.otu.tax$OTU,
                            "] ",
                            hm.otu.tax$Taxonomy,
                            sep='')

# Order samples according to Treatment and Compartment and OTUs by taxonomy
sample.order <- factor(env_16S$Compartment, levels = c("BS","AS","P")):factor(env_16S$Treatment, levels = c("Ctrl","TBAC","TEC","DBP"))
hm <- hm.tmp[order(sample.order),
             row.names(hm.otu.tax[order(hm.otu.tax$Taxonomy),])]

# Set palette and draw heatmap
mypal <- colorRampPalette(c("White", "Black"))

# svg('Heatmap_16S.svg')
heatmap(as.matrix(t(sqrt(hm))), revC = TRUE, 
        Colv = NA, Rowv = NA, scale = "none", col = mypal(256),
        labRow = hm.otu.tax$full.id)
# dev.off()

#######
# ITS #
#######

### Biomass - qPCR

# ANOVA with posthoc 
summary(aov(log2(qPCR$ITS.copies.g.material.1in100) ~ Plasticiser * Compartment, data = qPCR))
cld(lsmeans(
  aov(log2(qPCR$ITS.copies.g.material.1in100) ~ Treat.order * Comp.order, data = env_ITS), 
  ~ Treat.order * Comp.order), 
  Letters = letters, 
  adjust = "Tukey")

# Plot data from different treatment combinations as bars with 95% CI and individual data points
nTreatments <- 4
qPCR$i_Comp <- as.numeric(factor(qPCR$Compartment , levels = c('BS','AS','P')))
qPCR$i_Plasticiser <- as.numeric(factor(qPCR$Plasticiser , levels = c('Ctrl' , 'TBAC' , 'TEC' , 'DBP')))
qPCR$x_barplot <- (qPCR$i_Comp - 1) * (nTreatments + 1) + qPCR$i_Plasticiser + 0.5

bargraph.CI(env_ITS$Comp.order, log2(qPCR$ITS.copies.g.material.1in100), env_ITS$Treat.order,
            legend = FALSE, ylim = c(22,36), xlab = "", ylab = "ITS copy number per g material",
            main = "ITS copy number per g material", err.width = 0.05, col=plas.cols)
points(qPCR$x_barplot , log2(qPCR$ITS.copies.g.material.1in100), pch = 16 , cex = 0.4 , col = 'black')

### Alpha diversity

# ANOVA with Tukey HSD 
for(i in colnames(env_ITS[,5:7])){
  print(i)
  print(summary(aov(env_ITS[,i] ~ Treatment * Compartment, data = env_ITS)))
}

for(i in colnames(env_ITS[,5:7])) {
  print(i)
  print(cld(lsmeans(
    aov(env_ITS[,i] ~ Treatment * Compartment, data = env_ITS), 
    ~ Compartment * Treatment), 
    Letters = letters), 
    adjust = "Tukey")
}

nTreatments <- 4
env_ITS$i_Comp.order <- as.numeric(factor(env_ITS$Comp.order, 
                                          levels = c('a.BS','b.AS','c.P')))
env_ITS$i_Treatment <- as.numeric(factor(env_ITS$Treatment, 
                                         levels = c('Ctrl','TBAC','TEC','DBP')))
env_ITS$x_barplot <- (env_ITS$i_Comp.order - 1) * (nTreatments + 1) + env_ITS$i_Treatment + 0.5

for(i in colnames(env_ITS[,5:7])){
  bargraph.CI(env_ITS$Comp.order, env_ITS[,i],env_ITS$Treat.order,
              legend = FALSE, xlab = "", ylab = i, main = i, err.width = 0.05, 
              col = plas.cols , ylim = c(0 , 1.4 * max(env_ITS[,i])))
  points(env_ITS$x_barplot , env_ITS[,i] , pch = 16 , cex = 0.4 , col = 'black')
}

### Beta diversity

# PERMANOVA of Hellinger transformed OTUs
adonis2(sqrt(otu_ITS) ~ Treatment * Compartment, data = env_ITS, method = 'euc', permu = 4999)

### Visualise compositional differences with RDA

# Generate an RDA object constrained by Plasticiser:Compartment
otu_ITS.rda <- rda(sqrt(otu_ITS) ~ TreatComp, data = env_ITS)
ord = otu_ITS.rda

# Set colors, size and transparency for site symbols in RDA
addTrans <- function(col, alpha = 1) {
  rgb.mat <- col2rgb(col) / 255
  apply(rgb.mat, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha))
}

compartment.cex <- c(2, 4.5, 8)
cex = compartment.cex[factor(env_ITS$Comp.order)]

axis.percent <- function(ordination){
  round((100*eigenvals(ordination)[1:2]/ordination$tot.chi[[1]]),digits=2)
}

# Plot the RDA
scaling.val = 3
plot(ord, 
     type='n', scaling=scaling.val, 
     xlab=paste("RDA1 (",axis.percent(
       ord)[[1]],"%)",sep=""),
     ylab=paste("RDA2 (",axis.percent(
       ord)[[2]],"%)",sep=""))

points(ord, dis='sp', pch=4, col='grey', cex=0.6, scaling=scaling.val)

points(ord, dis='sites', pch=21, scaling = scaling.val,
       bg = addTrans(plas.cols[factor(env_ITS$Treat.order)], 180/255),
       cex = compartment.cex[factor(env_ITS$Comp.order)])

sd.val = 3
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1])),],
       pch=4,col='darkred',cex=0.6)	
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1]))),],
       pch=4,col='darkred',cex=0.6) 
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2])),],
       pch=4,col='darkred',cex=0.6)
points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2]))),],
       pch=4,col='darkred',cex=0.6) 
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1])),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1]))),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2])),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)
orditorp(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2]))),], 
         "sp", pch="+", col="black", air = 0.8, scaling = scaling.val)

legend("topleft",legend=unique(factor(env_ITS$Treatment)),pch=19,
       col=plas.cols[unique(factor(env_ITS$Treat.order))])

# Alternative plot for supplementary figures with different shapes for each compartment:
points(ord, dis='sites', 
       pch = c(21, 22, 24)[factor(env_ITS$Comp.order)],  # Different shapes
       scaling = scaling.val,
       bg = addTrans(plas.cols[factor(env_ITS$Treat.order)], 180/255),
       cex = 1.5)  # Keep size consistent

##############################
# Generate a heatmap for ITS #
##############################

# Calculate the mean relative abundance of each OTUs in each treatment combination
otu_ITS_means.tmp <- aggregate(otu_ITS, by = list(env_ITS$TreatComp), mean)
row.names(otu_ITS_means.tmp) <- otu_ITS_means.tmp[,1]
otu_ITS_means.tmp <- otu_ITS_means.tmp[,-1]

# Subset the OTU table to populations present at >=X% relative abundance
otus_for_heatmap <- colnames(otu_ITS_means.tmp[,which(apply(otu_ITS_means.tmp,2,max) >= 0.01)])
hm.tmp <- otu_ITS[, otus_for_heatmap] 

# Get the taxonomy
hm.otu.tax <- taxonomy_ITS[colnames(hm.tmp),]

hm.otu.tax$full.id <- paste("[",
                            hm.otu.tax$OTU,
                            "] ",
                            hm.otu.tax$Taxonomy,
                            sep='')

# Order samples according to Treatment and Compartment and OTUs by taxonomy
sample.order <- factor(env_ITS$Compartment, levels = c("BS","AS","P")):factor(env_ITS$Treatment, levels = c("Ctrl","TBAC","TEC","DBP"))
hm <- hm.tmp[order(sample.order),
             row.names(hm.otu.tax[order(hm.otu.tax$Taxonomy),])]
hm <- hm[,!(colnames(hm) == "Otu38")] # OTU38 was passed by ITSx as a fungal sequence but is not recognised in Unite10; NCBI nr indicates it may be a protist

# Set palette and draw heatmap
mypal <- colorRampPalette(c("White", "Black"))

# svg('Heatmap_ITS.svg')
heatmap(as.matrix(t(sqrt(hm))), revC = TRUE, 
        Colv = NA, Rowv = NA, scale = "none", col = mypal(256),
        labRow = hm.otu.tax$full.id)
# dev.off()

