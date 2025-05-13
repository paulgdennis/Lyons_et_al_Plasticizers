### SCRIPT 1. Lyons et al., (2025) 

# Paul Dennis (p.dennis@uq.edu.au) and Rebecca Lyons (r.lyons@uq.edu.au)

### List of scripts ###
# Script 1: Data input and sub-setting
# Script 2: Physicochemical analyses
# Script 3: Biomass, alpha, beta diversity
# Script 4: Indicator analyses
# Script 5: Differential abundance testing

##############################
# 1 Data input and sub-setting

################
# Get 16S data #
################

# Metadata including alpha diversity
env_16S <- read.table("Data/env_16S_5000.csv", header = TRUE, sep = ',', row.names = 1)

# Set factor levels
env_16S$Treatment <- factor(env_16S$Treatment, levels = c("Ctrl", "TBAC", "TEC","DBP"))
env_16S$Compartment <- factor(env_16S$Compartment, levels = c("BS", "AS", "P"))

# Rarefied OTU table with gg2 taxonomy
otu_16S_all.tmp <- read.table('Data/otu_16S_with_GG2_tax_5000.csv', header = TRUE, sep = ',', row.names = 1)
otu_16S <- as.data.frame(t(otu_16S_all.tmp[,-61])/5000)

# Check the samples in the OTU table and metadata are in the same order
all(row.names(env_16S) == row.names(otu_16S))

# Store the taxonomy in a new object
taxonomy_16S <- data.frame(OTU = row.names(otu_16S_all.tmp), 
                           Taxonomy = otu_16S_all.tmp[,61]) # Make taxonomy file
row.names(taxonomy_16S) <- taxonomy_16S$OTU

# Non-rarefied OTU table
otu_16S.nr <- read.table("Data/otu_16S_non_rarefied.csv", header = TRUE, sep = ',', row.names = 1)
all(rownames(env_16S) == colnames(otu_16S.nr))

# Weighted UniFrac
w_unifrac_16S.tmp <- read.table('Data/w_unifrac_16S_5000.csv', header = TRUE, sep = ',', row.names = 1)
w_unifrac_16S <- as.matrix(w_unifrac_16S.tmp)

# Check the samples in the weighted unifrac matrix and metadata are in the same order
all(row.names(w_unifrac_16S) == row.names(otu_16S))
all(colnames(w_unifrac_16S) == row.names(otu_16S))

#################
# Get ITS2 data #
#################

# Metadata including alpha diversity
env_ITS <- read.table("Data/env_ITS2_5850.csv", header = TRUE, sep = ',', row.names = 1)

# Set factor levels
env_ITS$Treatment <- factor(env_ITS$Treatment, levels = c("Ctrl", "TBAC", "TEC","DBP"))
env_ITS$Compartment <- factor(env_ITS$Compartment, levels = c("BS", "AS", "P"))

# OTU table
otu_ITS_all.tmp <- read.table('Data/otu_ITS2_with_Unite10_tax_5850.csv', header = TRUE, sep = ',', row.names = 1)
otu_ITS <- as.data.frame(t(otu_ITS_all.tmp[,-61])/5850)

# Check the samples in the OTU table and metadata are in the same order
all(row.names(env_ITS) == row.names(otu_ITS))

# Store the taxonomy in a new object
taxonomy_ITS <- data.frame(OTU = row.names(otu_ITS_all.tmp), 
                           Taxonomy = otu_ITS_all.tmp[,61]) # Make taxonomy file
row.names(taxonomy_ITS) <- taxonomy_ITS$OTU

# Non-rarefied OTU table
otu_ITS.nr <- read.table("Data/otu_ITS2_non_rarefied.csv", header = TRUE, sep = ',', row.names = 1)
all(rownames(env_ITS) == colnames(otu_ITS.nr))

#####################
# Get chemical data #
#####################

# GPC
chem_gpc <- read.table("Data/GPC.csv", header=TRUE, sep=',', row.names = 1)
chem_gpc$Plasticiser <- factor(chem_gpc$Plasticiser,
                               levels = c("Ctrl","TBAC","TEC","DBP"))

# NMR
chem_nmr <- read.table("Data/NMR.csv", header=TRUE, sep=',', row.names = 1)
chem_nmr$Plasticiser <- factor(chem_nmr$Plasticiser,
                               levels = c("Ctrl","TBAC","TEC","DBP"))

#################
# Get qPCR data #
#################

qPCR <- read.table('Data/qPCR.csv', header = TRUE, sep=',', row.names = 1)
qPCR$Plasticiser <- factor(qPCR$Plasticiser,
                           levels = c("Ctrl","TBAC","TEC","DBP"))

######################
# Colors and symbols #
######################

plas.cols <- c('#339e2bff','#b0de8aff','#a6cce3ff','#377EB8')

