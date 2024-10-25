# Load required libraries
library(Seurat)
library(monocle)
library(monocle3)
library(tidyverse)
library(patchwork)
library(viridis)
library(ggpubr)

# Clean Seurat object by removing 'Undefined' cluster
scRNA_cd8_clean <- scRNA_cd8[, scRNA_cd8$cluster != "Undefined"]

# Function to create CellDataSet (CDS) and preprocess it
create_cds <- function(seurat_obj) {
  data <- as(as.matrix(seurat_obj@assays$RNA@counts), 'sparseMatrix')
  pd <- new("AnnotatedDataFrame", data = seurat_obj@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  
  cds <- newCellDataSet(
    data,
    phenoData = pd,
    featureData = fd,
    expressionFamily = negbinomial.size()
  )
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds, cores = 8, relative_expr = FALSE)
  
  return(cds)
}

# Create and preprocess the CDS object
mycds <- create_cds(scRNA_cd8_clean)

# Filter high variable genes for trajectory analysis
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.075 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)

# Reduce dimensions using DDRTree
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree', norm_method = "none")

# Order cells along the trajectory
mycds <- orderCells(mycds)

# Define color palette for clusters
cluster_colors <- c("#F2AD00", "#FF0000", "#00A08A", "#F98400")

# Plot Trajectory - Cluster colored
plot_trajectory_cluster <- function(cds) {
  plot_cell_trajectory(cds, x = 1, y = 2, color_by = "cluster") +
    theme(
      legend.position = 'none',
      panel.border = element_blank()
    ) +
    scale_color_manual(values = cluster_colors)
}

# Plot Trajectory
plot_complex_trajectory <- function(cds) {
  plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "cluster") +
    scale_color_manual(values = cluster_colors) +
    theme(legend.title = element_blank())
}

# Save PDF for complex trajectory
pdf("Results/Fig4S_Pseudotimehierarchy.pdf", width = 4, height = 4, useDingbats = FALSE)
plot_complex_cell_trajectory(mycds, x = 1, y = 2, color_by = "cluster", show_branch_points = FALSE) +
  scale_color_manual(values = cluster_colors) +
  theme(legend.title = element_blank())
dev.off()

# Plot Pseudotime and Cluster Trajectories side-by-side
plot_pseudotime_trajectory <- function(cds) {
  p1 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "Pseudotime") +
    viridis::scale_colour_viridis(direction = 1) +
    theme(panel.border = element_blank())
  
  p2 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "cluster") +
    theme(panel.border = element_blank()) +
    scale_color_manual(values = cluster_colors)
  
  ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, align = "v")
}

# Save PDF of Pseudotime and Cluster Trajectories
pdf("Results/Fig3E_Pseudotime.pdf", width = 9, height = 4, useDingbats = FALSE)
print(plot_pseudotime_trajectory(mycds))
dev.off()

# Plot density distribution and trajectory colored by cluster
plot_density_and_trajectory <- function(cds) {
  data.peudo <- plot_cell_trajectory(cds)$data
  data.peudo$Group <- factor(data.peudo$Group, c("WT", "Tg"))
  
  p1 <- ggplot(data.peudo, aes(x = data_dim_1, color = Group, fill = Group)) +
    geom_density() +
    facet_wrap(vars(Group)) +
    theme(
      legend.position = 'none',
      panel.border = element_blank()
    ) +
    scale_color_manual(values = c('gray30', '#ec1c24')) +
    scale_fill_manual(values = c('gray70', '#FF9999')) +
    theme_classic()
  
  cds$Group <- factor(cds$Group, c("WT", "Tg"))
  
  p2 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "cluster") +
    theme(panel.border = element_blank()) +
    scale_color_manual(values = cluster_colors) +
    facet_wrap(vars(Group))
  
  ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, align = "h")
}

# Save PDF of Pseudotime Density Plot and Cluster Trajectory
pdf("Results/Fig3f_PseudotimeDensity.pdf", width = 8, height = 6, useDingbats = FALSE)
print(plot_density_and_trajectory(mycds))
dev.off()

# Extract data from trajectory
data.peudo <- plot_cell_trajectory(mycds)$data

# Set cluster levels
data.peudo$cluster <- factor(data.peudo$cluster, 
                             levels = c("NaiveCD8", "Early-activatedCD8", "EffectorCD8", "ExhaustedCD8"))

# Updated density plot with geom_smooth for better visualization
p1 <- ggplot(data.peudo, aes(x = data_dim_1, color = cluster, fill = cluster)) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +  # Using geom_smooth to better represent trends
  facet_wrap(vars(cluster)) +
  theme(legend.position = 'none', panel.border = element_blank()) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  theme_classic()

# Trajectory plot with updated colors
p2 <- plot_cell_trajectory(mycds, x = 1, y = 2, color_by = "cluster") +
  theme(legend.position = 'none', panel.border = element_blank()) +
  scale_color_manual(values = cluster_colors) +
  facet_wrap(vars(cluster))

# Save updated density and trajectory plots as PDF
pdf("Results/Fig4_PseudotimeDensity.pdf", width = 12, height = 8, useDingbats = FALSE)
ggarrange(p1, p2, ncol = 1, nrow = 2, align = "v")
dev.off()

# Changed geom_density to geom_smooth for better pseudotime visualization

pdf("Results/FigS6b_Pseudotime.pdf", width = 10, height = 12)

p_wave <- ggplot(data.peudo, aes(x = Pseudotime, color = cluster, fill = cluster)) +
  geom_smooth(method = "loess", se = TRUE, size = 1, alpha = 0.3) +  # Smooth trend for wave visualization
  facet_wrap(vars(cluster), ncol = 1) +
  theme(legend.position = 'none', panel.border = element_blank()) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  theme_classic()

print(p_wave)
dev.off()
