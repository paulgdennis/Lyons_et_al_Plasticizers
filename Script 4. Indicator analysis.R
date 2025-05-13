### SCRIPT 3. Lyons et al., (2025) 

# Paul Dennis (p.dennis@uq.edu.au) and Rebecca Lyons (r.lyons@uq.edu.au)

### List of scripts ###
# Script 1: Data input and sub-setting
# Script 2: Physicochemical analyses
# Script 3: Biomass, alpha, beta diversity
# Script 4: Indicator analyses
# Script 5: Differential abundance testing

# Indicator analysis

source('Script 1. Data input and subsetting.R')
library(labdsv)
library(ggplot2)
library(reshape2)
library(dplyr)

###################################
# Identify indicator taxa

#######
# 16S #
#######

### Run indval 
indval_res <- indval(x = sqrt(otu_16S[,which(apply(otu_16S*100, 2, max) >= 0.5)]), 
                     clustering = factor(env_16S$TreatComp))

# Get the group with the highest IndVal for each species
max_indval <- apply(indval_res$indval, 1, max)

# Filter by significance
sig <- indval_res$pval <= 0.05
sig_indvals <- max_indval[sig]

### Plot the top 10 most significant OTUs >=0.5% relative abundance for each treatment combination

# Create summary dataframe
indval_df <- as.data.frame(indval_res$indval)
indval_df$pval <- indval_res$pval
indval_df$OTU <- rownames(indval_df)

# Get list of treatment combinations (groups)
groups <- colnames(indval_df)[!colnames(indval_df) %in% c("pval", "OTU")]

# Create empty list to store plots
plot_list <- list()

for (grp in groups) {
  # Top 10 indicator OTUs for this group, filtering by p-value
  top_otus <- indval_df %>%
    filter(pval <= 0.01) %>%
    arrange(desc(.data[[grp]])) %>%
    slice_head(n = 10) %>%
    pull(OTU)
  
  # Subset OTU table and metadata
  df <- otu_16S[, top_otus, drop = FALSE]
  df$TreatComp <- env_16S$TreatComp
  
  # Keep only samples from this group
  df <- df[df$TreatComp == grp, ]
  
  # Long format
  df_melt <- melt(df, id.vars = "TreatComp")
  
  # Plot for this group
  p <- ggplot(df_melt, aes(x = variable, y = value*100)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = paste("Top 10 Indicator OTUs for", grp),
         x = "", y = "Relative abundance (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save to list
  plot_list[[grp]] <- p
}

# Show the plots one-by-one
plot_list[[1]]
plot_list[[2]]
plot_list[[3]]
plot_list[[4]]
plot_list[[5]]
plot_list[[6]]
plot_list[[7]]
plot_list[[8]]
plot_list[[9]]
plot_list[[10]]
plot_list[[11]]
plot_list[[12]]

# Create a directory to store plots and generate svg format files
#dir.create("indicator_svg_plots_16S", showWarnings = FALSE)

# #for (grp in names(plot_list)) {
# ggsave(filename = paste0("indicator_svg_plots_16S/Top10_Indicator_OTUs_", grp, ".svg"),
#        plot = plot_list[[grp]],
#        width = 8, height = 5,
#        device = "svg")
# }

### Generate a non-redundant list of the OTUs in the plots with their taxonomy
taxonomy_16S$OTU <- as.character(taxonomy_16S$OTU)

taxonomy_16S$Label <- sapply(strsplit(taxonomy_16S$Taxonomy, ";"), function(x) {
  label <- tail(x[!is.na(x) & x != ""], 1)
  ifelse(is.null(label) || label == "", "Unclassified", label)
})

# Combine top OTUs from each group into one list
top_otus_all <- lapply(groups, function(grp) {
  indval_df %>%
    filter(pval <= 0.05) %>%
    arrange(desc(.data[[grp]])) %>%
    slice_head(n = 10) %>%
    pull(OTU)
})

# Flatten and deduplicate
unique_otus <- unique(unlist(top_otus_all))

# Create a data frame with taxonomy
otu_tax_table <- taxonomy_16S %>%
  filter(OTU %in% unique_otus) %>%
  select(OTU, Taxonomy, Label)

# View or write to file
print(otu_tax_table)
#write.csv(otu_tax_table, "Top_Indicator_OTUs_Taxonomy.csv", row.names = FALSE)


#######
# ITS #
#######

### Run indval 
indval_res <- indval(x = sqrt(otu_ITS[,which(apply(otu_ITS*100, 2, max) >= 0.5)]), 
                     clustering = factor(env_ITS$TreatComp))

# Get the group with the highest IndVal for each species
max_indval <- apply(indval_res$indval, 1, max)

# Filter by significance
sig <- indval_res$pval <= 0.05
sig_indvals <- max_indval[sig]

### Plot the top 10 most significant OTUs >=0.5% relative abundance for each treatment combination

# Create summary dataframe
indval_df <- as.data.frame(indval_res$indval)
indval_df$pval <- indval_res$pval
indval_df$OTU <- rownames(indval_df)

# Get list of treatment combinations (groups)
groups <- colnames(indval_df)[!colnames(indval_df) %in% c("pval", "OTU")]

# Create empty list to store plots
plot_list <- list()

for (grp in groups) {
  # Top 10 indicator OTUs for this group, filtering by p-value
  top_otus <- indval_df %>%
    filter(pval <= 0.01) %>%
    arrange(desc(.data[[grp]])) %>%
    slice_head(n = 10) %>%
    pull(OTU)
  
  # Subset OTU table and metadata
  df <- otu_ITS[, top_otus, drop = FALSE]
  df$TreatComp <- env_ITS$TreatComp
  
  # Keep only samples from this group
  df <- df[df$TreatComp == grp, ]
  
  # Long format
  df_melt <- melt(df, id.vars = "TreatComp")
  
  # Plot for this group
  p <- ggplot(df_melt, aes(x = variable, y = value*100)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = paste("Top 10 Indicator OTUs for", grp),
         x = "", y = "Relative abundance (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save to list
  plot_list[[grp]] <- p
}

# Show the plots one-by-one
plot_list[[1]]
plot_list[[2]]
plot_list[[3]]
plot_list[[4]]
plot_list[[5]]
plot_list[[6]]
plot_list[[7]]
plot_list[[8]]
plot_list[[9]]
plot_list[[10]]
plot_list[[11]]
plot_list[[12]]

# # Create a directory to store plots and generate svg format files
# dir.create("indicator_svg_plots_ITS2", showWarnings = FALSE)
# 
# #for (grp in names(plot_list)) {
# ggsave(filename = paste0("indicator_svg_plots_ITS2/Top10_Indicator_OTUs_", grp, ".svg"),
#        plot = plot_list[[grp]],
#        width = 8, height = 5,
#        device = "svg")
# }

### Generate a non-redundant list of the OTUs in the plots with their taxonomy
taxonomy_ITS$OTU <- as.character(taxonomy_ITS$OTU)

taxonomy_ITS$Label <- sapply(strsplit(taxonomy_ITS$Taxonomy, ";"), function(x) {
  label <- tail(x[!is.na(x) & x != ""], 1)
  ifelse(is.null(label) || label == "", "Unclassified", label)
})

# Combine top OTUs from each group into one list
top_otus_all <- lapply(groups, function(grp) {
  indval_df %>%
    filter(pval <= 0.05) %>%
    arrange(desc(.data[[grp]])) %>%
    slice_head(n = 10) %>%
    pull(OTU)
})

# Flatten and deduplicate
unique_otus <- unique(unlist(top_otus_all))

# Create a data frame with taxonomy
otu_tax_table <- taxonomy_ITS %>%
  filter(OTU %in% unique_otus) %>%
  select(OTU, Taxonomy, Label)

# View or write to file
print(otu_tax_table)
#write.csv(otu_tax_table, "Top_Indicator_OTUs_ITS_Taxonomy.csv", row.names = FALSE)

