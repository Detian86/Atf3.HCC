# Load necessary libraries
library(Seurat)
library(GSVA)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# Extract subset of pathways of interest from GSVA results
selected_pathways <- c(
  "KEGG_TRYPTOPHAN_METABOLISM",
  "KEGG_PYRUVATE_METABOLISM",
  "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_N_GLYCAN_BIOSYNTHESIS",
  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  "GOBP_FATTY_ACID_BETA_OXIDATION",
  "GOBP_FATTY_ACID_BIOSYNTHETIC_PROCESS",
  "GOBP_FATTY_ACID_ELONGATION",
  "KEGG_CITRATE_CYCLE_TCA_CYCLE",
  "GOBP_ARGININE_METABOLIC_PROCESS",
  "KEGG_ARGININE_AND_PROLINE_METABOLISM"
)

# Subset GSVA results
gsva_sub <- gsva[selected_pathways, ]
gsva_sub <- t(gsva_sub)

# Check column names consistency between scRNA and GSVA results
if (!all(colnames(scRNA_macro) == rownames(gsva_sub))) {
  stop("Column names of scRNA_macro do not match with row names of gsva_sub.")
}

# Add pathway scores to Seurat metadata
scRNA_macro@meta.data <- cbind(scRNA_macro@meta.data, gsva_sub)
colnames(scRNA_macro@meta.data)[(ncol(scRNA_macro@meta.data) - ncol(gsva_sub) + 1) : ncol(scRNA_macro@meta.data)] <- colnames(gsva_sub)

# Dot plot visualization of pathway scores across clusters and groups
dotplot <- DotPlot(scRNA_macro, features = colnames(gsva_sub), group.by = "cluster", split.by = "Group", cluster.idents = FALSE) +
  theme(axis.text.x = element_text(angle = 90))

# Extract the data from DotPlot for further transformation
dotplot_data <- dotplot$data

# Transform the dotplot data into matrices for further visualization
# Create matrix for average expression
exp_mat <- dcast(dotplot_data[, c("features.plot", "id", "avg.exp")], features.plot ~ id)
row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[, -1] %>% as.matrix()

# Get quantiles for the expression matrix to understand data distribution
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

# Define color scheme for the heatmap
color_palette <- rev(brewer.pal(11, "Spectral"))[c(2, 4, 6, 8, 11)]
col_fun <- circlize::colorRamp2(c(-0.28, -0.08, 0, 0.02, 0.11), c("#3288BD", "#ABDDA4", "#FFFFBF", "#FDAE61", "#E31A1C"))

# Generate the heatmap using ComplexHeatmap
Heatmap(
  exp_mat,
  name = "Pathway Score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  row_names_gp = gpar(fontsize = 5),
  border = "black",
  column_title = "Clustered Dotplot"
)
