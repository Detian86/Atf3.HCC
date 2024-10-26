# Load necessary libraries
library(Seurat)
library(monocle)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(ggplot2)

# GSVA Analysis
# Define the expression matrix
expr <- as.matrix(scRNA_cd8@assays$RNA@data)

# Check that column names match between expression matrix and pseudotime metadata
if (!exists("mycds") || !("Pseudotime" %in% colnames(pData(mycds)))) {
  stop("Pseudotime metadata not found in 'mycds'. Ensure 'mycds' is defined and contains pseudotime information.")
}

if (!all(colnames(expr) == names(pData(mycds)$Pseudotime))) {
  stop("Column names do not match between expression matrix and pseudotime metadata.")
}

# Define pathway gene sets from MSigDB
msgdH <- msigdbr(species = "Mus musculus", category = "H")
msgdC2 <- msigdbr(species = "Mus musculus", category = "C2")
msgdC5 <- msigdbr(species = "Mus musculus", category = "C5")

# Convert pathway data to gene sets
geneSetH <- split(msgdH$gene_symbol, msgdH$gs_name)
geneSetC2 <- split(msgdC2$gene_symbol, msgdC2$gs_name)
geneSetC5 <- split(msgdC5$gene_symbol, msgdC5$gs_name)

# Filter pathways of interest from C2 and C5 collections based on selected keywords
selected_keywords <- c("GLYCOLYSIS", "OXIDATIVE_PHOSPHORYLATION", "FATTY", "TCA", "PENTOSE", "ARGININE", "PURINE")

geneSetC2 <- geneSetC2[grepl(paste(selected_keywords, collapse = "|"), names(geneSetC2))]
geneSetC5 <- geneSetC5[grepl(paste(selected_keywords, collapse = "|"), names(geneSetC5))]

# Perform GSVA for each gene set collection
gsvaH <- gsva(expr, gset.idx.list = geneSetH, kcdf = "Gaussian", method = "gsva", parallel.sz = 1)
gsvaC2 <- gsva(expr, gset.idx.list = geneSetC2, kcdf = "Gaussian", method = "gsva", parallel.sz = 1)
gsvaC5 <- gsva(expr, gset.idx.list = geneSetC5, kcdf = "Gaussian", method = "gsva", parallel.sz = 1)

# Combine GSVA results with pseudotime data
# Convert GSVA results into data frames for merging with pseudotime
gsvaH_df <- data.frame(t(gsvaH))
gsvaH_df$sample_name <- rownames(gsvaH_df)

gsvaC2_df <- data.frame(t(gsvaC2))
gsvaC2_df$sample_name <- rownames(gsvaC2_df)

gsvaC5_df <- data.frame(t(gsvaC5))
gsvaC5_df$sample_name <- rownames(gsvaC5_df)

# Merge GSVA data frames with pseudotime data
data.pseudo.gsvaH <- merge(data.peudo, gsvaH_df, by = "sample_name")
data.pseudo.gsvaC2 <- merge(data.peudo, gsvaC2_df, by = "sample_name")
data.pseudo.gsvaC5 <- merge(data.peudo, gsvaC5_df, by = "sample_name")

# Function to plot GSVA results over pseudotime for a given pathway
plot_pathway_activity <- function(data, pathway) {
  if (!pathway %in% colnames(data)) {
    stop(paste("Pathway", pathway, "not found in the provided data."))
  }
  
  ggplot(data, aes(x = Pseudotime, y = .data[[pathway]], color = Group, fill = Group)) +
    geom_smooth(method = "loess") +
    scale_color_manual(values = c('#ec1c24', 'black')) +
    scale_fill_manual(values = c('#FF9999', 'gray70')) +
    theme_bw() +
    theme(
      axis.title.x = element_text(colour = "black"),
      axis.title.y = element_text(colour = "black"),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black")
    )
}

# Example: Save the plot for a specific pathway from C2 collection
pathway_of_interest <- "WP_AEROBIC_GLYCOLYSIS"

if (pathway_of_interest %in% colnames(data.pseudo.gsvaC2)) {
  pdf_filename <- paste0("Results/", pathway_of_interest, "_over_pseudotime.pdf")
  pdf(pdf_filename, width = 12, height = 8, useDingbats = FALSE)
  print(plot_pathway_activity(data.pseudo.gsvaC2, pathway_of_interest))
  dev.off()
} else {
  message(paste("Pathway", pathway_of_interest, "not found in the data. Skipping plot generation."))
}
