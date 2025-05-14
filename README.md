# Lyons_et_al_Plasticizers
This repository contains the code used to analyse the 16S and ITS2 amplicon data reported in our associated manuscript: 

Rebecca Lyons, Clement Matthew Chan, Andrew R. Parry, Catherine M. E. Hodal, Jiarui Sun, Paul Lant, Steven Pratt, Bronwyn Laycock, Paul G. Dennis (2025) Relative to a common phthalate, citrate-based plasticizers exert minimal impact on plastisphere bacterial community composition during biopolymer biodegradation. Environmental Science and Technology. Under revision.

*Corresponding author: p.dennis@uq.edu.au

ABSTRACT 
Citrate-based plasticizers are often regarded as being less environmentally harmful alternatives to phthalates like dibutyl phthalate (DBP). Here, we tested these assertions by investigating the response of soil microbial communities to biopolymer samples plasticized with acetyl tributyl citrate (ATBC), triethyl citrate (TEC), DBP, or nothing (control). Samples were buried in soil for 14 weeks, after which the biomass and diversity of bacterial and fungal communities associated with the biopolymer, attached soil and bulk soil were characterized using metabarcoding and quantitative PCR. Differences between buried (incubated) and non-incubated biopolymer samples were also analyzed using X-ray micro-computed tomography, gel permeation chromatography and 1H-nuclear magnetic resonance spectroscopy to assess biopolymer degradation and plasticizer migration. Surface degradation and plasticizer migration were observed for all treatments, with minor impacts of plasticizers on microbial biomass and alpha diversity. Critically, however, plastisphere bacterial communities formed two groups that were compositionally distinct. One group comprised the DBP samples while the other comprised the control, ATBC and TEC samples. This indicates that DBP can impact soil bacterial community composition on release from biopolymer blends, with potential consequences for ecosystem function, while the effects of the citrate-based plasticizers were reduced, supporting their use as less environmentally impactful alternatives to phthalates. 

The repository has the following analytical R scripts:
- Script 1: Data input and sub-setting
- Script 2: Physicochemical analyses
- Script 3: Biomass, alpha, beta diversity
- Script 4: Indicator analyses
- Script 5: Differential abundance testing

In addition, the repository includes the following data files:
- GPC.csv 
- NMR.csv 
- qPCR.csv 
- env_16S_5000.csv
- otu_16S_with_GG2_tax_5000.csv
- otu_16S_non_rarefied.csv
- env_ITS2_5850.csv
- otu_ITS2_with_Unite10_tax_5850.csv
- otu_ITS2_non_rarefied.csv

Lastly, there is a text file with commands for our bioinformatics pipeline.
