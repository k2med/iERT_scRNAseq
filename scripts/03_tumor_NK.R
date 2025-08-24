# ---------------------------
# Title: Fig.3 — NK subclustering, UMAP, marker dotplots, OR heatmap; Fig.4A; Fig.S3/S4
# Purpose:
#   - Re-cluster NK cells (Harmony), visualize fine NK subsets (UMAP)
#   - Dotplots (NK markers; NK function markers) with source data
#   - Odds-ratio heatmap (subset enrichment per sample) with source data
#   - Polar barplots (cytokines per group) with per-group source data
#   - GSVA (hallmark/GO/KEGG/REACTOME) per NK subset heatmap with source data
# Inputs:
#   - ../data/Tumor_scRNA-seq_seurat.rds
#   - custom_colors.R  (defines nk_umap_colors, etc.)
#   - plot_heatmap.R   (heatmap helpers: heatmap_simple, top_markers_left, rectangle_annotation_coordinates)
#   - Fig4A_IREA_source_data.txt (for Fig.4A)
# Outputs:
#   - Fig.3A.umap.pdf
#   - Fig.3A.dotplot.pdf                 + Fig3A_dotplot_source_data.txt
#   - Fig.S3C.pdf                        + FigS3C_dotplot_source_data.txt
#   - Fig.S3A.pdf                        (NK markers heatmap)
#   - Fig.3B.pdf                         + Fig3B_OR_heatmap_source_data.txt
#   - Fig.4A/Fig.4A.<Group>.pdf          + Fig.4A/Fig4A_<Group>_source_data.txt
#   - Fig.S4A.pdf                        + FigS4A_source_data.txt
# ---------------------------

# Packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)
library(GSVA)
library(purrr)
library(egg)
library(grid)

source("custom_colors.R")
source("plot_heatmap.R")

# ---------------------------
# 0) Load dataset and subset NK cells
# ---------------------------
seu_all <- readRDS("../data/Tumor_scRNA-seq_seurat.rds")
seu_nk  <- subset(seu_all, subset = coarse_cell_type == "NK cell")

# ---------------------------
# 1) Build neighbor graph, cluster, and UMAP on Harmony dimensions
# ---------------------------
dims_use    <- 1:50
res_cluster <- 0.3

seu_nk <- FindNeighbors(seu_nk, reduction = "harmony", dims = dims_use)
seu_nk <- FindClusters(seu_nk, resolution = res_cluster)
seu_nk <- RunUMAP(seu_nk, reduction = "harmony", dims = dims_use, min.dist = 0.3)

# ---------------------------
# 2) Annotate clusters as fine NK cell types
# ---------------------------
seu_nk@meta.data <- seu_nk@meta.data |>
  mutate(
    fine_cell_type = case_when(
      seurat_clusters %in% c(0) ~ "Ccl5+ NK",
      seurat_clusters %in% c(1) ~ "Ccr2+ NK",
      seurat_clusters %in% c(2) ~ "Gzmg+ NK",
      seurat_clusters %in% c(3) ~ "Xcl1+ NK",
      seurat_clusters %in% c(4) ~ "Mki67+ NK"
    )
  )
seu_nk$fine_cell_type <- factor(
  seu_nk$fine_cell_type,
  levels = c("Xcl1+ NK","Ccl5+ NK","Ccr2+ NK","Gzmg+ NK","Mki67+ NK")
)
Idents(seu_nk) <- seu_nk$fine_cell_type

# ---------------------------
# 3) UMAP colored by fine NK cell types (Fig.3A.umap)
# ---------------------------
plot_nk_umap <- DimPlot(
  seu_nk,
  group.by = "fine_cell_type",
  pt.size  = 1,
  cols     = nk_umap_colors,
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

ggsave("Fig.3A.umap.pdf", plot_nk_umap, width = 24, height = 20, units = "cm")

umap_df <- data.frame(cell_id = colnames(seu_nk),
                      seu_nk@reductions$umap@cell.embeddings,
                      fine_cell_type = seu_nk$fine_cell_type
)

write.table(umap_df, "Fig3A_umap_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# 4) Dotplot of NK marker genes (Fig.3A.dotplot)
# ---------------------------
markers <- c("Ncr1","Fcgr3","Xcl1","Rora","Ccl5","Ccl4","Ccr2","Klrb1a","Gzmg","Gzmf","Mki67","Top2a")
gene <- intersect(markers, rownames(seu_nk))

plot_dot_raw <- DotPlot(seu_nk, features = gene)
dot_df <- plot_dot_raw$data[, c("id","features.plot","pct.exp","avg.exp.scaled")]

write.table(dot_df, "Fig3A_dotplot_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_nk_dot <- ggplot(dot_df, aes(x = id, y = features.plot)) +
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
  "Fig.3A.dotplot.pdf",
  egg::set_panel_size(
    plot_nk_dot,
    width  = grid::unit(length(gene) / 2.5, "cm"),
    height = grid::unit(length(levels(seu_nk$fine_cell_type)) / 2.5, "cm")
  ),
  width = 10, height = 10, units = "cm"
)

# ---------------------------
# 5) Dotplot of NK function genes (Fig.S3C)
# ---------------------------
markers <- c('Prf1','Gzma','Gzmb','Gzmc','Eomes','Tbx21','Rora','Gata3','Rorc','Runx3','Klf2','Tcf7',
             'Itga1','Itga2','Itgae','Itgam','Cxcr3','Cxcr4','Cxcr6','Cx3cr1','Ccr7','Ifng','Xcl1','Il5',
             'Il13','Ccl3','Ccl4','Ccl5','Nkg7','Klrc1','Klrk1','Klrg1','Cd96','Cd226','Lag3','Pdcd1',
             'Tigit','Mki67','Birc5')
gene <- intersect(markers, rownames(seu_nk))

plot_dot_raw <- DotPlot(seu_nk, features = gene)
dot_df <- plot_dot_raw$data[, c("id","features.plot","pct.exp","avg.exp.scaled")]

write.table(dot_df, "FigS3C_dotplot_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_nk_dot <- ggplot(dot_df, aes(x = id, y = features.plot)) +
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
  "Fig.S3C.pdf",
  egg::set_panel_size(
    plot_nk_dot,
    width  = grid::unit(length(gene) / 2.5, "cm"),
    height = grid::unit(length(levels(seu_nk$fine_cell_type)) / 2.5, "cm")
  ),
  width = 20, height = 10, units = "cm"
)

# ---------------------------
# 6) NK marker heatmap (Fig.S3A)
# Steps: MAST → significant → per-gene top log2FC → order → shuffle cells → draw heatmap
# ---------------------------
markers_all <- FindAllMarkers(
  seu_nk,
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

meta_nk <- seu_nk@meta.data
meta_nk_shuffled <- meta_nk %>%
  dplyr::mutate(row_id = rownames(meta_nk)) %>%
  dplyr::group_by(fine_cell_type) %>%
  dplyr::group_modify(~ dplyr::sample_n(.x, nrow(.x))) %>%
  dplyr::ungroup() %>%
  as.data.frame()
rownames(meta_nk_shuffled) <- meta_nk_shuffled$row_id

expr_mat <- as.matrix(
  seu_nk@assays$RNA@data[rownames(markers_sig_unique_ordered),
                         rownames(meta_nk_shuffled)]
)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Cell_type = meta_nk_shuffled$fine_cell_type,
  col = list(Cell_type = nk_umap_colors),
  annotation_name_gp = grid::gpar(fontface = "plain")
)

plot_markers_heatmap <- heatmap_simple(
  expr_mat,
  top_annotation = top_annotation,
  column_title   = "NK cells",
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

pdf("Fig.S3A.pdf", width = 10, height = 12)
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
  meta_nk_shuffled$fine_cell_type
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
# 7) Odds-ratio heatmap (Fig.3B) + source data
# ---------------------------
sample_levels <- c(
  "Ctrl_D04_Center","LDRT_D04_Center","SBRT_D04_Center","ERT_D04_Center",
  "Ctrl_D04_Peri","LDRT_D04_Peri","SBRT_D04_Peri","ERT_D04_Peri"
)
fine_nk_levels <- c("Xcl1+ NK","Ccl5+ NK","Gzmg+ NK","Ccr2+ NK","Mki67+ NK")

odds_ratio_mat <- matrix(
  nrow = length(fine_nk_levels),
  ncol = length(sample_levels),
  dimnames = list(fine_nk_levels, sample_levels)
)
pval_mat <- odds_ratio_mat

meta_all <- seu_all@meta.data
meta_all$fine_cell_type <- "Other"
meta_all[rownames(seu_nk@meta.data), "fine_cell_type"] <- as.character(seu_nk$fine_cell_type)

for (fine_type in rownames(odds_ratio_mat)) {
  for (sample_id in colnames(odds_ratio_mat)) {
    contingency_table <- as.data.frame.array(
      table(meta_all$fine_cell_type == fine_type,
            meta_all$orig.ident      == sample_id)
    )
    fisher_res <- fisher.test(contingency_table)
    odds_ratio_mat[fine_type, sample_id] <- fisher_res$estimate
    pval_mat[fine_type,  sample_id]      <- fisher_res$p.value
  }
}

pval_adj_mat <- matrix(
  p.adjust(pval_mat, method = "BH"),
  nrow = nrow(pval_mat), ncol = ncol(pval_mat),
  dimnames = dimnames(pval_mat)
)

odds_ratio_mat <- odds_ratio_mat[fine_nk_levels, sample_levels, drop = FALSE]
pval_adj_mat   <- pval_adj_mat[fine_nk_levels, sample_levels, drop = FALSE]

sig_symbols <- ifelse(pval_adj_mat < 1e-4, "****",
                      ifelse(pval_adj_mat < 1e-3, "***",
                             ifelse(pval_adj_mat < 1e-2, "**",
                                    ifelse(pval_adj_mat < 0.05, "*", ""))))

col_fun   <- circlize::colorRamp2(c(0, 1, 3), c("#62bde6","white","#de7379"))
col_split <- rep(seq_len(ncol(odds_ratio_mat) / 4), each = 4, length.out = ncol(odds_ratio_mat))

pdf("Fig.3B.pdf", width = 6, height = 10)
ComplexHeatmap::Heatmap(
  odds_ratio_mat,
  col = col_fun,
  name = "Odds ratio",
  cluster_columns = FALSE,
  cluster_rows    = FALSE,
  show_row_dend   = FALSE,
  rect_gp         = grid::gpar(col = "#dddddd"),
  show_heatmap_legend = TRUE,
  width  = grid::unit(8 * 0.3 + 0.1 * 2, "cm"),
  height = grid::unit(5 * 0.3, "cm"),
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp    = grid::gpar(fontsize = 6),
  column_split    = col_split,
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(sig_symbols[i, j], x = x, y = y, gp = grid::gpar(fontsize = 4))
  }
)
dev.off()

# Source data (long format)
or_df <- as.data.frame(odds_ratio_mat)
or_df$cell <- rownames(or_df)
or_long <- pivot_longer(or_df, cols = -cell, names_to = "sample", values_to = "OR")

padj_df <- as.data.frame(pval_adj_mat)
padj_df$cell <- rownames(padj_df)
padj_long <- pivot_longer(padj_df, cols = -cell, names_to = "sample", values_to = "padj")

or_source_data <- left_join(or_long, padj_long, by = c("cell", "sample"))
write.table(or_source_data, "Fig3B_OR_heatmap_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# 8) Fig.4A — Polar barplots of cytokines per group (+ per-group source data)
# ---------------------------
input_df <- read.table("Fig4A_IREA_source_data.txt", header = TRUE, sep = "", check.names = FALSE)

prep_one_group <- function(df) {
  df %>%
    arrange(padj) %>%
    slice_head(n = 20) %>%
    arrange(ES) %>%
    mutate(
      Cytokine   = factor(Cytokine, levels = Cytokine),
      angle      = 90 - 360 / n() * (row_number() - 0.5),
      angle      = if_else(angle < -90, angle + 180, angle),
      log10padj  = -log10(padj) * sign(ES),
      log10padj  = pmax(pmin(log10padj, 5), -5)   # clip to [-5, 5]
    )
}

group_list <- split(input_df, input_df$Group)
prep_list  <- imap(group_list, ~ prep_one_group(.x) %>% mutate(Group = .y))
df_top     <- bind_rows(prep_list)

es_min <- min(df_top$ES, na.rm = TRUE)
es_max <- max(df_top$ES, na.rm = TRUE)

plot_polar <- function(df_one_group) {
  ggplot(df_one_group, aes(x = Cytokine)) +
    geom_col(aes(y = ES, fill = log10padj)) +
    scale_fill_gradientn(
      colors = c("#62bde6","white","#de7379"),
      limits = c(-5, 5),
      name   = "-log10(adjusted p value)"
    ) +
    geom_text(aes(y = es_max + 600, label = Cytokine, angle = angle),
              size = 6 / .pt) +
    coord_polar() +
    ylim(c(es_min - 200, es_max + 800)) +
    theme_void() +
    theme(
      text = element_text(size = 6),
      legend.position = "bottom",
      legend.title    = element_text(size = 6),
      legend.text     = element_text(size = 6)
    )
}

dir.create("Fig.4A", showWarnings = FALSE)

walk(prep_list, function(d) {
  g <- plot_polar(d)
  ggsave(
    filename = file.path("Fig.4A", paste0("Fig.4A.", unique(d$Group), ".pdf")),
    egg::set_panel_size(g, width = grid::unit(5, "cm"), height = grid::unit(5, "cm")),
    width = 10, height = 15, units = "cm"
  )
  write.table(
    d %>% select(Group, Cytokine, ES, padj, log10padj),
    file = file.path("Fig.4A", paste0("Fig4A_", unique(d$Group), "_source_data.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
})

# ---------------------------
# 9) Fig.S4A — GSVA heatmap of NK subsets (hallmark/GO/KEGG/REACTOME)
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

expr_matrix <- as.matrix(seu_nk@assays$RNA@data)

gsva_results <- gsva(
  expr          = expr_matrix,
  gset.idx.list = genesets,
  kcdf          = "Gaussian",
  method        = "gsva"
)

gsva_by_cell <- as.data.frame(t(gsva_results))   # cells x pathways
gsva_by_cell$cell_type <- seu_nk$fine_cell_type

pathways_of_interest <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "GOBP_NATURAL_KILLER_CELL_ACTIVATION",
  "GOBP_NATURAL_KILLER_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_NATURAL_KILLER_CELL_CYTOKINE_PRODUCTION",
  "GOBP_NATURAL_KILLER_CELL_DEGRANULATION",
  "GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL",
  "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
  "GOBP_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION"
)

cluster_pathway_mat <- gsva_by_cell %>%
  dplyr::select(all_of(pathways_of_interest), cell_type) %>%
  group_by(cell_type) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

rownames(cluster_pathway_mat) <- cluster_pathway_mat$cell_type
cluster_pathway_mat <- as.matrix(t(cluster_pathway_mat[, -1]))  # pathways x clusters

cluster_pathway_mat_z <- t(scale(t(cluster_pathway_mat)))       # row-wise Z

pdf("Fig.S4A.pdf", width = 8, height = 6)
Heatmap(
  cluster_pathway_mat_z,
  col = c("#62bde6","white","#de7379"),
  name = "Average score",
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  rect_gp = grid::gpar(col = "#dddddd"),
  show_heatmap_legend = TRUE,
  width  = grid::unit(5 * 0.3, "cm"),
  height = grid::unit(9 * 0.3, "cm"),
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp    = grid::gpar(fontsize = 6)
)
dev.off()

write.table(cluster_pathway_mat_z, "FigS4A_source_data.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)
