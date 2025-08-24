# ---------------------------
# Title: Fig.5 — cDC subclustering, UMAP, marker dotplot, heatmap, and GSVA violin plot
# Purpose:
#   - Re-cluster cDC cells using Harmony dimensions
#   - Visualize fine cDC subsets by UMAP
#   - Plot cDC marker dotplot (with source data)
#   - Draw cDC marker heatmap (Fig.S5A)
#   - Compute GSVA and draw violin plot for selected pathways (with source data)
# Inputs:
#   - ../data/Tumor_scRNA-seq_seurat.rds
#   - custom_colors.R         (defines cdc_umap_colors, heatmap_colors, etc.)
#   - plot_heatmap.R          (heatmap_simple, top_markers_left, rectangle_annotation_coordinates)
# Outputs:
#   - Fig.5A.umap.pdf
#   - Fig.5A.dotplot.pdf                  + Fig5A_dotplot_source_data.txt
#   - Fig.S5A.pdf                         (cDC marker heatmap)
#   - Fig.5B.pdf                           + Fig5B_vlnplot_source_data.txt
# ---------------------------

# Packages
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)
library(GSVA)
library(egg)
library(grid)

source("custom_colors.R")
source("plot_heatmap.R")

# ---------------------------
# 0) Load dataset and subset cDC cells
# ---------------------------
seu_all <- readRDS("../data/Tumor_scRNA-seq_seurat.rds")
seu_cdc  <- subset(seu_all, subset = coarse_cell_type == "cDC")

# ---------------------------
# 1) Build neighbor graph, cluster, and UMAP on Harmony dimensions
# ---------------------------
dims_use    <- 1:50
res_cluster <- 0.5

seu_cdc <- FindNeighbors(seu_cdc, reduction = "harmony", dims = dims_use)
seu_cdc <- FindClusters(seu_cdc, resolution = res_cluster)
seu_cdc <- RunUMAP(seu_cdc, reduction = "harmony", dims = dims_use, min.dist = 0.3)

# ---------------------------
# 2) Annotate clusters as fine cDC cell types
# ---------------------------
seu_cdc@meta.data <- seu_cdc@meta.data |>
  mutate(
    fine_cell_type = case_when(
      (seurat_clusters %in% c(0)) ~ "Cd209a+ cDC2",
      (seurat_clusters %in% c(1)) ~ "Xcr1+ cDC1",
      (seurat_clusters %in% c(2)) ~ "Ccr7+ cDC1",
      (seurat_clusters %in% c(3)) ~ "Mki67+ cDC"
    )
  )
seu_cdc$fine_cell_type <- factor(
  seu_cdc$fine_cell_type,
  levels = c("Xcr1+ cDC1","Ccr7+ cDC1","Cd209a+ cDC2","Mki67+ cDC")
)
Idents(seu_cdc) <- seu_cdc$fine_cell_type

# ---------------------------
# 3) UMAP colored by fine cDC cell types (Fig.5A.umap)
# ---------------------------
plot_cdc_umap <- DimPlot(
  seu_cdc,
  group.by = "fine_cell_type",
  pt.size  = 1,
  cols     = cdc_umap_colors,
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

ggsave("Fig.5A.umap.pdf", plot_cdc_umap, width = 24, height = 20, units = "cm")

umap_df <- data.frame(cell_id = colnames(seu_cdc),
                      seu_cdc@reductions$umap@cell.embeddings,
                      fine_cell_type = seu_cdc$fine_cell_type
)

write.table(umap_df, "Fig5A_umap_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# 4) Dotplot of cDC marker genes (Fig.5A.dotplot)
# ---------------------------
markers <- c("Xcr1","Clec9a","Itgae","Ccr7","Ccl5","Ccl22","Cd209a","Cd300a","Ccr5","Mki67","Top2a","Stmn1")
gene <- intersect(markers, rownames(seu_cdc))

plot_dot_raw <- DotPlot(seu_cdc, features = gene)
dot_df <- plot_dot_raw$data[, c("id","features.plot","pct.exp","avg.exp.scaled")]

write.table(dot_df, "Fig5A_dotplot_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_cdc_dot <- ggplot(dot_df, aes(x = id, y = features.plot)) +
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
  "Fig.5A.dotplot.pdf",
  egg::set_panel_size(
    plot_cdc_dot,
    width  = grid::unit(length(gene) / 2.5, "cm"),
    height = grid::unit(length(levels(seu_cdc$fine_cell_type)) / 2.5, "cm")
  ),
  width = 10, height = 10, units = "cm"
)

# ---------------------------
# 5) cDC marker heatmap (Fig.S5A)
# Steps: MAST → significant → per-gene top log2FC → order → shuffle cells → draw heatmap
# ---------------------------
markers_all <- FindAllMarkers(
  seu_cdc,
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

meta_cdc <- seu_cdc@meta.data
meta_cdc_shuffled <- meta_cdc %>%
  dplyr::mutate(row_id = rownames(meta_cdc)) %>%
  dplyr::group_by(fine_cell_type) %>%
  dplyr::group_modify(~ dplyr::sample_n(.x, nrow(.x))) %>%
  dplyr::ungroup() %>%
  as.data.frame()
rownames(meta_cdc_shuffled) <- meta_cdc_shuffled$row_id

expr_mat <- as.matrix(
  seu_cdc@assays$RNA@data[rownames(markers_sig_unique_ordered),
                         rownames(meta_cdc_shuffled)]
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Cell_type = meta_cdc_shuffled$fine_cell_type,
  col = list(Cell_type = cdc_umap_colors),
  annotation_name_gp = grid::gpar(fontface = "plain")
)

plot_markers_heatmap <- heatmap_simple(
  expr_mat,
  top_annotation = top_annotation,
  column_title   = "cDC cells",
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

pdf("Fig.S5A.pdf", width = 10, height = 12)
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
  meta_cdc_shuffled$fine_cell_type
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
# 6) Fig.5B — GSVA vlnplot of cDC subsets
# ---------------------------
get_genesets <- function(species = "Mus musculus", category, subcategory = NULL) {
  msigdbr(species = species, category = category, subcategory = subcategory) %>%
    dplyr::select(gs_name, gene_symbol) %>%
    split(x = .$gene_symbol, f = .$gs_name)
}

genesets <- c(
  get_genesets(category = "H"),
  get_genesets(category = "C5", subcategory = "GO:BP"),
  get_genesets(category = "C2", subcategory = "KEGG"),
  get_genesets(category = "C2", subcategory = "CP:REACTOME")
)

expr_matrix <- as.matrix(seu_cdc@assays$RNA@data)

gsva_results <- gsva(
  expr          = expr_matrix,
  gset.idx.list = genesets[pathways_of_interest],
  kcdf          = "Gaussian",
  method        = "gsva"
)

gsva_by_cell <- as.data.frame(t(gsva_results))   # cells x pathways
gsva_by_cell$cell_type <- seu_cdc$fine_cell_type

# keep POIs and convert to long
gsva_by_cell_sub <- gsva_by_cell[, c(pathways_of_interest, "cell_type")]
gsva_long <- reshape2::melt(
  gsva_by_cell_sub,
  id.vars = "cell_type",
  variable.name = "pathway",
  value.name   = "score"
)

# export source data
write.table(gsva_long[,c("cell_type", "score")],
            "Fig5B_vlnplot_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# violin plot
p_vln <- ggplot(
  subset(gsva_long, pathway %in% pathways_of_interest),
  aes(x = cell_type, y = score, fill = cell_type)
) +
  geom_violin(scale = "width", size = 0.116) +
  geom_boxplot(outlier.shape = NA, size = 0.116, width = 0.12, fill = "white") +
  xlab("") + ylab("") +
  scale_fill_manual(values = cdc_umap_colors) +
  theme(
    text = element_text(size = 6),
    panel.grid.major = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor = element_blank(),
    panel.background  = element_blank(),
    axis.line  = element_line(colour = "black", size = 0.116),
    axis.ticks = element_line(colour = "black", size = 0.116),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(colour = "black", size = 6),
    legend.title = element_text(size = 6)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  "Fig.5B.pdf",
  egg::set_panel_size(p_vln, width = grid::unit(2.5, "cm"), height = grid::unit(2, "cm")),
  width = 8, height = 8, units = "cm"
)

p_vln <- p_vln +
  ggpubr::stat_pwc(
    aes(group = cell_type),
    label.size = 6 * 25.4 / 72,
    tip.length = 0,
    size = 0.116,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.adj.signif"
  )

ggsave(
  "Fig.5B.sig.pdf",
  egg::set_panel_size(p_vln, width = grid::unit(2.5, "cm"), height = grid::unit(3.5, "cm")),
  width = 8, height = 8, units = "cm"
)
