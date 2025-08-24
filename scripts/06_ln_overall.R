# ---------------------------
# Title: Fig.7A — LN overall UMAP
# Purpose:
#   - Draw overall UMAP colored by coarse cell type
# Inputs:
#   - ../data/LN_scRNA-seq_seurat.rds
#   - custom_colors.R  (defines overall_umap_colors)
# Outputs:
#   - Fig.7A.overall.pdf                 (overall UMAP)
# ---------------------------

library(Seurat)
library(ggplot2)
source("custom_colors.R")

# 0) Load
seu_all <- readRDS("../data/LN_scRNA-seq_seurat.rds")

# 1) Overall UMAP by coarse_cell_type
plot_overall <- DimPlot(
  seu_all,
  group.by = "coarse_cell_type",
  pt.size  = 0.2,
  cols     = overall_umap_colors,
  raster   = FALSE
) +
  labs(x = "UMAP1", y = "UMAP2", title = NULL) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text   = element_blank(),
    aspect.ratio = 1
  )

ggsave("Fig.7A.overall.pdf", plot_overall, width = 24, height = 20, units = "cm")

umap_df <- data.frame(cell_id = colnames(seu_all),
                      seu_all@reductions$umap@cell.embeddings,
                      coarse_cell_type = seu_all$coarse_cell_type
)

write.table(umap_df, "Fig7A_umap_overall_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
