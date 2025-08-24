# ---------------------------
# Title: Fig.6 — CD8+ T subclustering, UMAP, marker dotplots, heatmap, AUCell density, and TCR diversity
# Purpose:
#   - Re-cluster CD8+ T cells (Harmony) and visualize fine subsets (UMAP)
#   - Plot CD8+ T marker/function dotplots with source data (Fig.6A, Fig.S6B)
#   - Draw CD8+ T marker heatmap (Fig.S6A)
#   - Score gene programs by AUCell and plot density on UMAP (Fig.6B)
#   - Summarize TCR diversity and visualize clonal sizes (Fig.S6C–D)
# Inputs:
#   - ../data/Tumor_scRNA-seq_seurat.rds
#   - custom_colors.R  (defines: cd8t_umap_colors, overall_umap_colors, heatmap_colors)
#   - plot_heatmap.R   (heatmap_simple, top_markers_left, rectangle_annotation_coordinates)
#   - TCR files: <base_dir>/<orig.ident>/TCR_filtered_contig_annotations.csv
# Outputs:
#   - Fig.6A.umap.pdf
#   - Fig.6A.dotplot.pdf                  + Fig6A_dotplot_source_data.txt
#   - Fig.S6B.dotplot.pdf                 + FigS6B_dotplot_source_data.txt
#   - Fig.S6A.pdf                         (CD8+ T marker heatmap)
#   - Fig.6B.<signature>.pdf              (AUCell density per signature)
#   - Fig.S6C.pdf                         + FigS6C_shannon_source_data.txt
#   - Fig.S6D.pdf                         (UMAP split by clonal size)
# ---------------------------

# Packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(AUCell)
library(RColorBrewer)
library(Nebulosa)
library(scRepertoire)
library(stringr)
library(egg)
library(grid)
library(reshape2)

source("custom_colors.R")
source("plot_heatmap.R")

# ---------------------------
# 0) Load dataset and subset CD8+ T cells
# ---------------------------
seu_all <- readRDS("../data/Tumor_scRNA-seq_seurat.rds")
seu_cd8t  <- subset(seu_all, subset = coarse_cell_type == "CD8+ T cell")

# ---------------------------
# 1) Build neighbor graph, cluster, and UMAP on Harmony dimensions
# ---------------------------
dims_use    <- 1:50
res_cluster <- 0.5

seu_cd8t <- FindNeighbors(seu_cd8t, reduction = "harmony", dims = dims_use)
seu_cd8t <- FindClusters(seu_cd8t, resolution = res_cluster)
seu_cd8t <- RunUMAP(seu_cd8t, reduction = "harmony", dims = dims_use, min.dist = 0.3)

# ---------------------------
# 2) Annotate clusters as fine CD8+ T cell types
# ---------------------------
seu_cd8t@meta.data <- seu_cd8t@meta.data |>
  mutate(
    fine_cell_type = case_when(
      (seurat_clusters %in% c(0)) ~ "Effector-like CD8+ T (Tef)",
      (seurat_clusters %in% c(1)) ~ "Proliferating CD8+ T",
      (seurat_clusters %in% c(2)) ~ "ISG+ CD8+ T",
      (seurat_clusters %in% c(3)) ~ "Exhausted CD8+ T (Tex)",
      (seurat_clusters %in% c(5)) ~ "Precursor of exhausted CD8+ T (Tpex)",
      (seurat_clusters %in% c(6)) ~ "Naive/memory CD8+ T"
    )
  )
seu_cd8t <- subset(seu_cd8t, cells = rownames(seu_cd8t@meta.data[!is.na(seu_cd8t@meta.data$fine_cell_type), ]))
seu_cd8t$fine_cell_type <- factor(
  seu_cd8t$fine_cell_type,
  levels = c("Naive/memory CD8+ T",
             "Precursor of exhausted CD8+ T (Tpex)",
             "Effector-like CD8+ T (Tef)",
             "Exhausted CD8+ T (Tex)",
             "ISG+ CD8+ T",
             "Proliferating CD8+ T")
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

ggsave("Fig.6A.umap.pdf", plot_cd8t_umap, width = 24, height = 20, units = "cm")

umap_df <- data.frame(cell_id = colnames(seu_cd8t),
                      seu_cd8t@reductions$umap@cell.embeddings,
                      fine_cell_type = seu_cd8t$fine_cell_type
)

write.table(umap_df, "Fig6A_umap_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# 4) Dotplots (Fig.6A, Fig.S6B) + source data
# ---------------------------
## 4.1 Fig.6A marker set
markers <- c("Cd3d","Cd8a","Tcf7","Lef1","Sell","Ccr7","Bach2","Slamf6","Cxcr6","Gzmd","Prf1",
             "Pdcd1","Havcr2","Tigit","Lag3","Cxcr3","Gzmk","Ifit3","Mki67","Top2a")
gene <- intersect(markers, rownames(seu_cd8t))

plot_dot_raw <- DotPlot(seu_cd8t, features = gene)
dot_df <- plot_dot_raw$data[, c("id","features.plot","pct.exp","avg.exp.scaled")]

write.table(dot_df, "Fig6A_dotplot_source_data.txt",
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
  "Fig.6A.dotplot.pdf",
  egg::set_panel_size(
    plot_cd8t_dot,
    width  = grid::unit(length(gene) / 2.5, "cm"),
    height = grid::unit(length(levels(seu_cd8t$fine_cell_type)) / 2.5, "cm")
  ),
  width = 16, height = 10, units = "cm"
)

## 4.2 Fig.S6B function set
markers <- c("Tcf7","Tox","Tbx21","Eomes","Bcl2","Id2","Lef1","Pdcd1","Slamf6",
             "Havcr2","Cd101","Cxcr3","Cxcr6","Cx3xr1","S1pr1","S1pr5","Il2",
             "Tnf","Ifng","Gzma","Gzmb","Prf1")
gene <- intersect(markers, rownames(seu_cd8t))

plot_dot_raw <- DotPlot(seu_cd8t, features = gene)
dot_df <- plot_dot_raw$data[, c("id","features.plot","pct.exp","avg.exp.scaled")]

write.table(dot_df, "FigS6B_dotplot_source_data.txt",
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
  "Fig.S6B.pdf",
  egg::set_panel_size(
    plot_cd8t_dot,
    width  = grid::unit(length(gene) / 2.5, "cm"),
    height = grid::unit(length(levels(seu_cd8t$fine_cell_type)) / 2.5, "cm")
  ),
  width = 16, height = 10, units = "cm"
)

# ---------------------------
# 5) CD8+ T marker heatmap (Fig.S6A)
# Steps: MAST → significant → per-gene top log2FC → order → shuffle cells → draw heatmap
# ---------------------------
markers_all <- FindAllMarkers(
  seu_cd8t,
  assay = "RNA",
  only.pos = FALSE,
  test.use = "MAST",
  min.pct = 0.1,
  logfc.threshold = 0.1
)

markers_sig <- subset(markers_all, (p_val_adj < 0.05) & (avg_log2FC > 0.25))

markers_sig_unique <- markers_sig %>%
  dplyr::group_by(gene) %>%
  dplyr::slice_max(avg_log2FC, n = 1) %>%
  dplyr::ungroup()

markers_sig_unique_ordered <- as.data.frame(
  markers_sig_unique[order(markers_sig_unique$cluster,
                           markers_sig_unique$avg_log2FC,
                           decreasing = c(FALSE, TRUE)), ]
)
rownames(markers_sig_unique_ordered) <- markers_sig_unique_ordered$gene

meta_cd8t <- seu_cd8t@meta.data
meta_cd8t_shuffled <- meta_cd8t %>%
  dplyr::mutate(row_id = rownames(meta_cd8t)) %>%
  dplyr::group_by(fine_cell_type) %>%
  dplyr::group_modify(~ dplyr::sample_n(.x, nrow(.x))) %>%
  dplyr::ungroup() %>%
  as.data.frame()
rownames(meta_cd8t_shuffled) <- meta_cd8t_shuffled$row_id

expr_mat <- as.matrix(
  seu_cd8t@assays$RNA@data[rownames(markers_sig_unique_ordered),
                         rownames(meta_cd8t_shuffled)]
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Cell_type = meta_cd8t_shuffled$fine_cell_type,
  col = list(Cell_type = cd8t_umap_colors),
  annotation_name_gp = grid::gpar(fontface = "plain")
)

plot_markers_heatmap <- heatmap_simple(
  expr_mat,
  top_annotation = top_annotation,
  column_title   = "CD8+ T cells",
  name           = "Markers_heatmap",
  scale_rows     = TRUE,
  width          = grid::unit(3, "in"),
  height         = grid::unit(6, "in"),
  raster_quality = 7,
  legend_name    = "Relative expression",
  color_palette  = heatmap_colors,
  show_annotation_name = FALSE,
  color_range    = 2 * seq(-1, 1, length.out = length(heatmap_colors) - 1)
)

pdf("Fig.S6A.pdf", width = 10, height = 12)
suppressWarnings({
  ComplexHeatmap::draw(
    top_markers_left(expr_mat, top_n = 10,
                     ord = markers_sig_unique_ordered, ord_col = "cluster") + plot_markers_heatmap,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legends = FALSE
  )
})

rect_coords <- rectangle_annotation_coordinates(
  markers_sig_unique_ordered$cluster,
  meta_cd8t_shuffled$fine_cell_type
)

ComplexHeatmap::decorate_heatmap_body("Markers_heatmap", {
  grid::grid.rect(
    x = grid::unit(rect_coords$x, "native"),
    y = grid::unit(rect_coords$y, "native"),
    width  = grid::unit(rect_coords$w, "native"),
    height = grid::unit(rect_coords$h, "native"),
    hjust = 0, vjust = 1,
    gp = grid::gpar(col = "white", lty = 1, lwd = 1)
  )
})
dev.off()

# ---------------------------
# 6) Fig.6B — AUCell programs: density on UMAP (Nebulosa)
# ---------------------------
expr_matrix <- as.matrix(seu_cd8t@assays$RNA@data)

# Build rankings
cells_rankings <- AUCell_buildRankings(expr_matrix, nCores = 3, plotStats = FALSE)

gene_sets <- list(
  stemness      = c("Tcf7","Cxcr5","Cd28","Gzmk","Ccr7","Il7r","Bcl6","Sell","Cd27"),
  cytotoxicity  = c("Prf1","Ifng","Nkg7","Gzmb","Gzma","Klrk1","Klrb1","Ctsw","Cst7"),
  naiveness     = c("Ccr7","Tcf7","Lef1","Sell","Il7r"),
  exhaustion    = c("Lag3","Ctla4","Pdcd1","Havcr2","Tox","Icos","Tnfrsf4","Tigit")
)

cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)  # AUC matrix: signatures × cells

# Add as a new assay "AUC" so Nebulosa/FeaturePlot can read it as features
auc_mat <- as.matrix(cells_AUC@assays@data@listData$AUC)
auc_assay <- CreateAssayObject(data = auc_mat)
seu_cd8t_auc <- seu_cd8t
seu_cd8t_auc[["AUC"]] <- auc_assay

# Plot density per program on UMAP
program_order <- c("stemness","naiveness","cytotoxicity","exhaustion")
for (pg in program_order) {
  DefaultAssay(seu_cd8t_auc) <- "AUC"
  p_den <- Nebulosa::plot_density(
    seu_cd8t_auc, features = pg, reduction = "umap", size = 0.5
  ) +
    scale_color_gradientn(colors = rev(colorRampPalette(brewer.pal(9, "RdBu"))(30))) +
    labs(x = "UMAP1", y = "UMAP2", title = NULL) +
    theme(
      axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
      panel.background = element_blank()
    )
  ggsave(
    filename = paste0("Fig.6B.", pg, ".pdf"),
    plot = egg::set_panel_size(p_den, width = grid::unit(4, "cm"), height = grid::unit(4, "cm")),
    width = 8, height = 8, units = "cm"
  )
}

# ---------------------------
# 7) Fig.S6C–D — TCR diversity (scRepertoire)
# ---------------------------
base_dir <- "../../data"
id_vec   <- levels(seu_cd8t$orig.ident)

get_barcode <- function(v) sapply(strsplit(v, "_", fixed = TRUE), `[`, 4)

# Read per-sample TCR and filter to barcodes present in object
tcr_list <- lapply(id_vec, function(id) {
  fp    <- file.path(base_dir, id, "TCR_filtered_contig_annotations.csv")
  tcr   <- read.csv(fp, stringsAsFactors = FALSE)
  cells <- rownames(seu_cd8t@meta.data[seu_cd8t$orig.ident == id, ])
  bc    <- get_barcode(cells)
  tcr[tcr$barcode %in% bc, ]
})
names(tcr_list) <- id_vec

samples <- sub("^([^_]+)_([^_]+)_.*$", "\\1_\\2", id_vec)
IDs     <- sub("^([^_]+)_([^_]+)_([^_]+)$", "\\3", id_vec)

tcr_merged <- combineTCR(tcr_list, samples = samples, ID = IDs)

seu_cd8t_tcr <- combineExpression(
  tcr_merged, seu_cd8t,
  cloneCall  = "strict",
  proportion = TRUE,
  cloneSize  = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1)
)

md <- seu_cd8t_tcr@meta.data %>% filter(!is.na(cloneSize))

# Shannon index per orig.ident
shannon_df <- md %>%
  count(orig.ident, CTstrict, name = "n") %>%
  group_by(orig.ident) %>%
  mutate(p = n / sum(n)) %>%
  summarise(Index = -sum(p * log(p)), .groups = "drop") %>%
  mutate(
    Group = str_split_fixed(orig.ident, "_", 3)[, 1],
    Site  = str_c(str_split_fixed(orig.ident, "_", 3)[, 2],
                  str_split_fixed(orig.ident, "_", 3)[, 3], sep = "_"),
    Group = factor(Group, levels = c("Ctrl","LDRT","SBRT","ERT"))
  )

write.table(shannon_df, "FigS6C_shannon_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

p_div <- ggplot(shannon_df, aes(x = Site, y = Index, group = Group)) +
  geom_path(linetype = "dashed", color="grey") +
  geom_point(size = 2, aes(color=Group), show.legend = TRUE) +
  theme_bw(base_size = 20) +
  labs(x = NULL, y = "Shannon diversity index") +
  theme(
    text = element_text(size = 6),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", size = 0.125),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size = 0.125),
    axis.ticks = element_line(colour = "black", size = 0.125),
    axis.text.x = element_text(colour = "black", angle = 60, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(colour = "black", size = 6),
    legend.position = "right",
    legend.title = element_text(size = 10)
  ) +
  scale_color_manual(values = overall_umap_colors)

ggsave(
  filename = "Fig.S6C.pdf",
  plot = egg::set_panel_size(p_div, width = grid::unit(4, "cm"), height = grid::unit(4, "cm")),
  width = 10, height = 10, units = "cm"
)

seu_cd8t_tcr_d10 <- subset(seu_cd8t_tcr, time == "D10")
plot_cd8t_tcr_umap <- DimPlot(
  seu_cd8t_tcr_d10,
  group.by = "cloneSize",
  pt.size  = 1,
  cols = c(`Rare (0 < X <= 1e-04)`='#c6dbef', `Small (1e-04 < X <= 0.001)`='#709bf5',
           `Medium (0.001 < X <= 0.01)`='#81c47e', `Large (0.01 < X <= 0.1)`='#e38c49',
           `Hyperexpanded (0.1 < X <= 1)`='#dc462d'),
  na.value = "#d9d9d9",
  split.by = "group",
  ncol = 2,
  raster = FALSE
) +
  labs(x = "UMAP1", y = "UMAP2", title = NULL) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1
  )
ggsave("Fig.S6D.pdf", plot_cd8t_tcr_umap, width = 24, height = 20, units = "cm")


umap_df <- data.frame(cell_id = colnames(seu_cd8t_tcr_d10),
                      seu_cd8t_tcr_d10@reductions$umap@cell.embeddings,
                      cloneSize = seu_cd8t_tcr_d10$cloneSize,
                      group = seu_cd8t_tcr_d10$group
)
head(umap_df)
write.table(umap_df, "FigS6D_umap_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
