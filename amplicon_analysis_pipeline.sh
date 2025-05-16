#!/bin/bash

# By: Jiarui Sun (email: jiarui.sun@uq.edu.au)

# Please note this script only record the key commands used in our amplicon analysis pipeline,
# which was developed by Adam Skarshewski(GitHub: https://github.com/Askars), and now maintained
# by Jiarui Sun. 

#===============================================================================#
## Step 1: Import sequences into QIIME, demultiplex, and visualise data summary
#===============================================================================#
# // Import sequences into QIIME:
#$ qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path manifest.csv \
#      --output-path demux.qza --input-format SingleEndFastqManifestPhred33

# // Create QIIME summaries/visualisations:
#$ qiime demux summarize --i-data demux.qza --o-visualization demux.qzv

# // Summarise collated data:
#$ qiime tools export --input-path demux.qzv --output-path demultiplexing_reports

# // Extract demultiplexed sequences
#$ qiime tools export --input-path demux.qza --output-path output_dir 
### Then concatenate all sequences into 'combined.fastq'


#===============================================================================#                
## Step	2.1: Processing ITS sequencing data
#===============================================================================# 
# // Run ITSx: 
### Firstly all sequences from fastq files were merged and converted into fasta file 'combined.fasta', then:
# $ itsx -i combined.fasta -o combined_itsx -t F  # Then concatenate all fasta files into 'ITSx_output.fasta'

# // Extract reference db for chimera checking:
#$ qiime tools export --input-path unite_8.3_2021_05_10_dynamic_s_seq.qza --output-path tmp/unite_ref

# // Check for chimeras:
#$ usearch -uchime2_ref ITSx_output.fasta -db tmp/unite_ref/dna-sequences.fasta \
       -notmatched no_chimeras.fasta -threads 15 -strand plus -mode balanced
### Then clean up header lines and save to 'cleaned.fasta'

# *** Use UPARSE pipeline for denovo OTU pick ***

# // Remove duplicates for rep set generation:
#$ usearch -fastx_uniques cleaned.fasta -fastaout uniques.fasta -sizeout -relabel Uniq

# // Generate representative sequence set:
#$ usearch -cluster_otus uniques.fasta -otus rep_set.fasta -relabel Otu

# // Generate OTU table:
#$ usearch -otutab cleaned.fasta -otus rep_set.fasta -otutabout non_normalised_otu.tsv


#===============================================================================#
## Step 2.2: Processing 16S or 18S sequencing data
#===============================================================================# 
# // Perform denoising and OTU clustering (UPARSE):
## Strip 16S primer: AAACTYAAAKGAATTGRCGG
$ usearch -fastx_truncate combined.fastq -stripleft 20 -fastqout stripped.fastq # 20 is the length of primer.

# // Truncate to 250 bp:
#$ usearch -fastx_truncate stripped.fastq -trunclen 250 -fastqout trimmed.fastq

# // Quality filtering reads:
#$ usearch -fastq_filter trimmed.fastq -fastqout filtered.fastq -fastq_maxee 1.0

# // Remove duplicates for rep set generation:
#$ usearch -fastx_uniques filtered.fastq -fastaout uniques.fasta -sizeout -relabel Uniq

# // Generate representative sequence set:
#$ usearch -cluster_otus uniques.fasta -otus rep_set.fasta -relabel Otu

# // Generate OTU table:
#$ usearch -otutab trimmed.fastq -otus rep_set.fasta -otutabout non_normalised_otu.tsv


#===============================================================================#
## Step 3: Taxonomy assignment
#===============================================================================#
# // Convert OTU table to BIOM format:
#$ biom convert -i non_normalised_otu.tsv -o otu.biom --to-hdf5

# // Import OTU table into QIIME:
#$ qiime tools import --type 'FeatureTable[Frequency]' --input-path otu.biom --output-path table.qza

# // Import representative sequences into QIIME:
#$ qiime tools import --type 'FeatureData[Sequence]' --input-path rep_set.fasta --output-path rep-seqs.qza

# // Perform taxonomic assignment:
#$ qiime feature-classifier classify-consensus-blast --i-query rep-seqs.qza \
#    --i-reference-reads silva_138_SSURef_Nr99_seq.qza --i-reference-taxonomy silva_138_SSURef_Nr99_tax.qza \
#    --o-classification taxonomy_blast.qza --o-search-results top-hits.qza
## Note: DB files for ITS are: 'unite_8.3_2021_05_10_dynamic_s_seq.qza' and 'unite_8.3_2021_05_10_dynamic_s_tax.qza'

# // Export taxonomy annotation results:
#$ qiime tools export --input-path taxonomy_blast.qza --output-path output/taxonomy
# There will be a 'taxonomy.tsv' file within the output/taxonomy/ directory.
# If there are mitochondrial or chloroplast OTUs, remove them.


#===============================================================================#
## Step 4: Phylogenetic tree                 
#===============================================================================#
# // Performing multiple sequence alignment
#$ qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza

# // Mask alignments:
#$ qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza

# // Build de-novo tree:
#$ qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza

# // Midpoint rooting tree:
#$ qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
