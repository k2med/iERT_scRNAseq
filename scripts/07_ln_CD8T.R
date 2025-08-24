# ---------------------------
# Title: Fig.7 — CD8+ T subclustering, UMAP, and marker dotplot
# Purpose:
#   - Re-cluster CD8+ T cells (Harmony) and visualize fine subsets by UMAP (Fig.7A)
#   - Plot CD8+ T marker/function dotplot with source data (Fig.7B)
# Inputs:
#   - ../data/LN_scRNA-seq_seurat.rds
#   - custom_colors.R  (defines cd8t_umap_colors)
# Outputs:
#   - Fig.7A.cd8t.pdf                  (UMAP of CD8+ T subsets)
#   - Fig.7B.pdf                       + Fig7B_dotplot_source_data.txt
# ---------------------------

# Packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(egg)
library(grid)
library(reshape2)
source("custom_colors.R")

# ---------------------------
# 0) Load dataset and subset CD8+ T cells
# ---------------------------
seu_all <- readRDS("../data/LN_scRNA-seq_seurat.rds")
seu_cd8t  <- subset(seu_all, subset = coarse_cell_type == "CD8+ T cell")

# ---------------------------
# 1) Build neighbor graph, cluster, and UMAP on Harmony dimensions
# ---------------------------
dims_use    <- 1:50
res_cluster <- 0.1

seu_cd8t <- FindNeighbors(seu_cd8t, reduction = "harmony", dims = dims_use)
seu_cd8t <- FindClusters(seu_cd8t, resolution = res_cluster)
seu_cd8t <- RunUMAP(seu_cd8t, reduction = "harmony", dims = dims_use, min.dist = 0.3)

# ---------------------------
# 2) Annotate clusters as fine CD8+ T cell types
# ---------------------------
seu_cd8t@meta.data <- seu_cd8t@meta.data |>
  mutate(
    fine_cell_type = case_when(
      (seurat_clusters %in% c(0,2)) ~ "Naive/memory CD8+ T",
      (seurat_clusters %in% c(1)) ~ "Precursor of exhausted CD8+ T (Tpex)"
    )
  )
seu_cd8t <- subset(seu_cd8t, cells = rownames(seu_cd8t@meta.data[!is.na(seu_cd8t@meta.data$fine_cell_type), ]))
seu_cd8t$fine_cell_type <- factor(
  seu_cd8t$fine_cell_type,
  levels = c("Precursor of exhausted CD8+ T (Tpex)",
             "Naive/memory CD8+ T")
)
Idents(seu_cd8t) <- seu_cd8t$fine_cell_type

# ---------------------------
# 3) UMAP colored by fine CD8+ T cell types (Fig.6A.umap)
# ---------------------------
plot_cd8t_umap <- DimPlot(
  seu_cd8t,
  group.by = "fine_cell_type",
  pt.size  = 1,
  cols     = cd8t_umap_colors,
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

ggsave("Fig.7A.cd8t.pdf", plot_cd8t_umap, width = 24, height = 20, units = "cm")

umap_df <- data.frame(cell_id = colnames(seu_cd8t),
                      seu_cd8t@reductions$umap@cell.embeddings,
                      fine_cell_type = seu_cd8t$fine_cell_type
)

write.table(umap_df, "Fig7A_umap_cd8t_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# 4) Dotplots (Fig.7B) + source data
# ---------------------------
markers <- c("Ccr7","Tcf7","Lef1","Sell","Slamf6","Cxcr3",
             "Xcl1","Gzmk","Pdcd1","Havcr2","Tigit","Lag3")
gene <- intersect(markers, rownames(seu_cd8t))

plot_dot_raw <- DotPlot(seu_cd8t, features = gene)
dot_df <- plot_dot_raw$data[, c("id","features.plot","pct.exp","avg.exp.scaled")]

write.table(dot_df, "Fig7B_dotplot_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_cd8t_dot <- ggplot(dot_df, aes(x = id, y = features.plot)) +
  geom_point(
    aes(fill = avg.exp.scaled, size = pct.exp),
    color = "black", shape = 21, stroke = 0.01
  ) +
  xlab("") + ylab("") +
  scale_fill_gradientn(colors = c("#62bde6","white","#de7379")) +
  scale_size(range = c(0, 3.2), limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +
  theme(
    text = element_text(size = 10),
    panel.grid.major = element_line(colour = "grey90", size = 0.25),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line  = element_line(colour = "black", size = 0.25),
    axis.ticks = element_line(colour = "black", size = 0.25),
    axis.text.x = element_text(colour = "black", angle = 60, vjust = 1, hjust = 1, size = 6, face = "italic"),
    axis.text.y = element_text(colour = "black", size = 6),
    legend.position = "bottom",
    legend.title    = element_text(size = 10)
  ) +
  guides(
    size = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 1, byrow = TRUE,
                        override.aes = list(stroke = 0.4)),
    fill = guide_colourbar(title.position = "top", title.hjust = 0.5)
  ) +
  coord_flip() +
  scale_x_discrete(limits = rev) 

ggsave(
  "Fig.7B.pdf",
  egg::set_panel_size(
    plot_cd8t_dot,
    width  = grid::unit(length(gene) / 2.5, "cm"),
    height = grid::unit(length(levels(seu_cd8t$fine_cell_type)) / 2.5, "cm")
  ),
  width = 16, height = 10, units = "cm"
)
