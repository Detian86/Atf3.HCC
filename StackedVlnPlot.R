# Load necessary libraries
library(Seurat)
library(MySeuratWrappers)
library(ggplot2)

# Define marker genes for plotting
Marker <- c("Ptprc", "Cd3g", "Cd3e", "Cd4", "Cd8a", "Cd8b1", "Cd19",
            "Ms4a1", "Mzb1", "Ncr1", "Klrb1c", "Tcrg-C1", "Itgam",
            "Cd14", "Fcgr3", "Itgax", "Flt3", "Siglech", "Lair1",
            "S100a8", "S100a9", "Top2a", "Mki67")

# Set cell identity levels in reverse order for better visualization
scRNA@active.ident <- factor(scRNA@active.ident,
                             levels = rev(c("Macrophages", "mDCs", "pDCs", "Neutrophil",
                                            "B", "Plasma", "CD4T", "CD8T", "NKT", "gdT", 
                                            "CyclingT", "NK", "Undetermined")))

# Define color palette
colors <- rev(c('#57C3F3', '#DCC1DD', '#476D87', '#58A4C3', '#E59CC4',
                         '#BD956A', '#D6E7A3', '#E95C59', '#F3B1A0', '#625D9E',
                         '#F1BB72', '#53A85F', '#E0D4CA'))
                         
# Define a function to generate the violin plot
create_violin_plot <- function(seurat_obj, features, colors) {
  MySeuratWrappers::VlnPlot(
    object = seurat_obj,
    features = features,
    direction = "horizontal",
    adjust = 1,
    cols = colors,
    stacked = TRUE,
    pt.size = 0
  ) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1))
}

# Create and display the stacked violin plot
violin_plot <- create_violin_plot(scRNA, Marker, colors)
print(violin_plot)