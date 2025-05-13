### SCRIPT 5. Lyons et al., (2025) 

# Paul Dennis (p.dennis@uq.edu.au) and Rebecca Lyons (r.lyons@uq.edu.au)

### List of scripts ###
# Script 1: Data input and sub-setting
# Script 2: Physicochemical analyses
# Script 3: Biomass, alpha, beta diversity
# Script 4: Indicator analyses
# Script 5: Differential abundance testing

# Differential abundance testing with DESeq2, Supporting Information

source('Script 1. Data input and subsetting.R')
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Set updated color mapping for treatments 
colors_treatments <- c(Ctrl = '#33a02c', TBAC = '#b2df8a', TEC = '#a6cee3', DBP = '#1f78b4')

# Define compartments and treatments
compartments <- c("BS", "AS", "P")
treatments <- c("Ctrl", "TBAC", "TEC", "DBP")

#######
# 16S #
#######

# Calculate mean relative abundance and keep OTU >= 0.5% relative abundance
otu_rel_abund <- sweep(otu_16S.nr, 2, colSums(otu_16S.nr), "/")
otu_mean_abund <- rowMeans(otu_rel_abund)
keep_otus <- names(otu_mean_abund[otu_mean_abund >= 0.005])

# Loop through each compartment, subset data, run DESeq2
dds_list <- list()
res_list <- list()

for (comp in compartments) {
  keep_samples <- env_16S$Compartment == comp
  meta_sub <- env_16S[keep_samples, , drop = FALSE]
  
  if (length(unique(meta_sub$Treatment)) < 2 || any(table(meta_sub$Treatment) < 2)) {
    message("[Skipping] Not enough replicates in compartment: ", comp)
    next
  }
  
  otu_sub <- otu_16S.nr[, keep_samples]
  
  dds <- DESeqDataSetFromMatrix(countData = otu_sub,
                                colData = meta_sub,
                                design = ~ Treatment)
  dds <- DESeq(dds)
  dds_list[[comp]] <- dds
  
  for (treat in treatments[treatments != "Ctrl"]) {
    res <- results(dds, contrast = c("Treatment", treat, "Ctrl"))
    res_df <- as.data.frame(res) %>%
      rownames_to_column(var = "OTU") %>%
      filter(OTU %in% keep_otus) %>%
      mutate(Treatment = treat,
             Compartment = comp,
             Significant = padj < 0.05 & !is.na(padj),
             Label = ifelse(Significant & OTU %in% keep_otus, OTU, NA))
    res_list[[paste(treat, comp, sep = ".")]] <- res_df
    
    # Export CSV
    outdir <- "deseq2_by_compartment_16S"
    dir.create(outdir, showWarnings = FALSE)
    write.csv(res_df, file = file.path(outdir, paste0("DESeq2_", treat, "_vs_Ctrl_in_", comp, ".csv")), row.names = FALSE)
  }
}

# Combine all results into one data frame for multi-treatment plots
all_res <- bind_rows(res_list)

# --- 3-panel volcano plots per compartment with all treatments ---
plots <- list()

for (comp in compartments) {
  df <- all_res %>% filter(Compartment == comp)
  
  if (nrow(df) == 0) {
    warning("No data available for volcano plot in compartment: ", comp)
    p <- ggplot() +
      theme_void() +
      ggtitle(paste("No data for", comp))
  } else {
    df$Treatment <- factor(df$Treatment, levels = names(colors_treatments))
    
    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Treatment)) +
      geom_point(alpha = 0.7) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_text_repel(aes(label = Label), size = 3, max.overlaps = 15) +
      scale_color_manual(values = colors_treatments) +
      theme_minimal(base_size = 14) +
      labs(title = paste("Volcano Plot:", comp), x = "Log2 Fold Change", y = "-Log10 adjusted p-value")
    
    ggsave(file.path("deseq2_by_compartment_16S", paste0("Volcano_All_Treatments_in_", comp, ".svg")), plot = p, width = 6, height = 5)
  }
  
  plots[[comp]] <- p
}

plots[["BS"]]
plots[["AS"]]
plots[["P"]]

#######
# ITS #
#######

# Calculate mean relative abundance and keep OTU >= 0.5% relative abundance
otu_rel_abund <- sweep(otu_ITS.nr, 2, colSums(otu_ITS.nr), "/")
otu_mean_abund <- rowMeans(otu_rel_abund)
keep_otus <- names(otu_mean_abund[otu_mean_abund >= 0.005])

# Loop through each compartment, subset data, run DESeq2
dds_list <- list()
res_list <- list()

for (comp in compartments) {
  keep_samples <- env_ITS$Compartment == comp
  meta_sub <- env_ITS[keep_samples, , drop = FALSE]
  
  if (length(unique(meta_sub$Treatment)) < 2 || any(table(meta_sub$Treatment) < 2)) {
    message("[Skipping] Not enough replicates in compartment: ", comp)
    next
  }
  
  otu_sub <- otu_ITS.nr[, keep_samples]
  
  dds <- DESeqDataSetFromMatrix(countData = otu_sub,
                                colData = meta_sub,
                                design = ~ Treatment)
  dds <- DESeq(dds)
  dds_list[[comp]] <- dds
  
  for (treat in treatments[treatments != "Ctrl"]) {
    res <- results(dds, contrast = c("Treatment", treat, "Ctrl"))
    res_df <- as.data.frame(res) %>%
      rownames_to_column(var = "OTU") %>%
      filter(OTU %in% keep_otus) %>%
      mutate(Treatment = treat,
             Compartment = comp,
             Significant = padj < 0.05 & !is.na(padj),
             Label = ifelse(Significant & OTU %in% keep_otus, OTU, NA))
    res_list[[paste(treat, comp, sep = ".")]] <- res_df
    
    # Export CSV
    outdir <- "deseq2_by_compartment_ITS2"
    dir.create(outdir, showWarnings = FALSE)
    write.csv(res_df, file = file.path(outdir, paste0("DESeq2_", treat, "_vs_Ctrl_in_", comp, ".csv")), row.names = FALSE)
  }
}

# Combine all results into one data frame for multi-treatment plots
all_res <- bind_rows(res_list)

# --- 3-panel volcano plots per compartment with all treatments ---
plots <- list()

for (comp in compartments) {
  df <- all_res %>% filter(Compartment == comp)
  
  if (nrow(df) == 0) {
    warning("No data available for volcano plot in compartment: ", comp)
    p <- ggplot() +
      theme_void() +
      ggtitle(paste("No data for", comp))
  } else {
    df$Treatment <- factor(df$Treatment, levels = names(colors_treatments))
    
    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Treatment)) +
      geom_point(alpha = 0.7) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_text_repel(aes(label = Label), size = 3, max.overlaps = 15) +
      scale_color_manual(values = colors_treatments) +
      theme_minimal(base_size = 14) +
      labs(title = paste("Volcano Plot:", comp), x = "Log2 Fold Change", y = "-Log10 adjusted p-value")
    
    ggsave(file.path("deseq2_by_compartment_ITS2", paste0("Volcano_All_Treatments_in_", comp, ".svg")), plot = p, width = 6, height = 5)
  }
  
  plots[[comp]] <- p
}

plots[["BS"]]
plots[["AS"]]
plots[["P"]]
