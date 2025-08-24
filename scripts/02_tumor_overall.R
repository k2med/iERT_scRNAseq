# ---------------------------
# Title: Fig.2 — Overall UMAP & site-faceted NK/CD8; Fig.S2C dotplot
# Purpose:
#   - Overall UMAP by coarse cell type
#   - Site-faceted UMAP for NK (top) and CD8 (bottom), colored by group
#   - Canonical marker dotplot with source data
# Inputs:
#   - ../data/Tumor_scRNA-seq_seurat.rds
#   - custom_colors.R  (defines overall_umap_colors)
# Outputs:
#   - Fig.2F.pdf                         (overall UMAP)
#   - Fig.2G.pdf                         (NK on top + CD8 at bottom, faceted by site)
#   - Fig.S2C.pdf                        + FigS2C_dotplot_source_data.txt
#   - Fig.5D.pdf                         + FigS5D_dotplot_source_data.txt
# ---------------------------


library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
source("custom_colors.R")

# 0) Load
seu_all <- readRDS("../data/Tumor_scRNA-seq_seurat.rds")

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

ggsave("Fig.2F.pdf", plot_overall, width = 24, height = 20, units = "cm")

# 2) Subset
seu_cd8 <- subset(seu_all, subset = coarse_cell_type == "CD8+ T cell")
seu_nk  <- subset(seu_all, subset = coarse_cell_type == "NK cell")

# Global UMAP axis range for consistent panels
umap_all <- Embeddings(seu_all, reduction = "umap")
x_range  <- range(umap_all[, 1], na.rm = TRUE)
y_range  <- range(umap_all[, 2], na.rm = TRUE)

# Helper: build meta (UMAP + site + group)
build_meta <- function(seu_obj) {
  emb <- Embeddings(seu_obj, reduction = "umap") |>
    as.data.frame() |>
    rownames_to_column("cell_id")
  md  <- seu_obj@meta.data |>
    as.data.frame() |>
    rownames_to_column("cell_id")
  out <- left_join(md, emb, by = "cell_id")
  out$group <- factor(out$group, levels = c("ERT", "SBRT", "LDRT", "Ctrl"))
  out$site <- paste(
    sapply(strsplit(as.character(out$orig.ident), "_"), `[`, 2),
    sapply(strsplit(as.character(out$orig.ident), "_"), `[`, 3),
    sep = "_"
  )
  out
}

meta_cd8 <- build_meta(seu_cd8)
meta_nk  <- build_meta(seu_nk)

# Helper: site-faceted scatter (color by group)
make_site_plot <- function(df) {
  ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = group), size = 0.2) +
    scale_color_manual(values = overall_umap_colors) +
    facet_grid(~ site) +
    labs(x = "UMAP1", y = "UMAP2", title = NULL) +
    theme_classic() +
    theme(
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      strip.background = element_rect(color = "white", fill = "white"),
      aspect.ratio     = 1
    ) +
    scale_x_continuous(limits = x_range) +
    scale_y_continuous(limits = y_range)
}

plot_nk_site  <- make_site_plot(meta_nk)
plot_cd8_site <- make_site_plot(meta_cd8)

# 3) NK (top) + CD8 (bottom)
plot_fig2g <- plot_nk_site / plot_cd8_site
ggsave("Fig.2G.pdf", plot_fig2g, width = 60, height = 40, units = "cm")

umap_df <- data.frame(cell_id = colnames(seu_all),
                      seu_all@reductions$umap@cell.embeddings,
                      coarse_cell_type = seu_all$coarse_cell_type,
                      group = seu_all$group,
                      site = paste(
                        sapply(strsplit(as.character(seu_all$orig.ident), "_"), `[`, 2),
                        sapply(strsplit(as.character(seu_all$orig.ident), "_"), `[`, 3),
                        sep = "_"
                      )
)

write.table(umap_df, "Fig2F_umap_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# 4) Dotplot of canonical marker genes (Fig.S2C)
# ---------------------------
Idents(seu_all) <- seu_all$coarse_cell_type

markers <- c("Vim","Krt20","Birc5",
             "Ptprc","Cd3d","Cd3e", "Cd8a", "Cd4", 
             "Ncr1","Klrd1","Tyrobp","Fcgr3",
             "Lyz2","Apoe","C1qa",
             "S100a9","S100a8","G0s2",
             "Cd74","H2-Aa","Ccr7",
             "Klk1","Siglech","Irf8",
             "Cd79a","Ms4a1","Cd19",
             "Col3a1", "Col1a1", 
             "Pecam1", "Cdh5", "Cd93")
gene <- intersect(markers, rownames(seu_all))

plot_dot_raw <- DotPlot(seu_all, features = gene)
dot_df <- plot_dot_raw$data[, c("id","features.plot","pct.exp","avg.exp.scaled")]

write.table(dot_df, "FigS2C_dotplot_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_all_dot <- ggplot(dot_df, aes(x = id, y = features.plot)) +
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
  "Fig.S2C.pdf",
  egg::set_panel_size(
    plot_all_dot,
    width  = grid::unit(length(gene) / 2.5, "cm"),
    height = grid::unit(length(levels(seu_all$coarse_cell_type)) / 2.5, "cm")
  ),
  width = 20, height = 16, units = "cm"
)

# ---------------------------
# 5) Dotplot of Xcl1 and Xcr1 (Fig.5D)
# ---------------------------
Idents(seu_all) <- seu_all$coarse_cell_type

markers <- c("Xcl1","Xcr1")
gene <- intersect(markers, rownames(seu_all))

plot_dot_raw <- DotPlot(seu_all, features = gene)
dot_df <- plot_dot_raw$data[, c("id","features.plot","pct.exp","avg.exp.scaled")]

write.table(dot_df, "Fig5D_dotplot_source_data.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_all_dot <- ggplot(dot_df, aes(x = id, y = features.plot)) +
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
    axis.text.y = element_text(colour = "black", size = 6, face = "italic"),
    axis.text.x = element_text(colour = "black", angle = 60, vjust = 1, hjust = 1, size = 6),
    legend.position = "bottom",
    legend.title    = element_text(size = 10)
  ) +
  guides(
    size = guide_legend(title.position = "top", title.hjust = 0.5, ncol = 1, byrow = TRUE,
                        override.aes = list(stroke = 0.4)),
    fill = guide_colourbar(title.position = "top", title.hjust = 0.5)
  ) +
  scale_y_discrete(limits = rev) 

ggsave(
  "Fig.5D.pdf",
  egg::set_panel_size(
    plot_all_dot,
    width  = grid::unit(length(levels(seu_all$coarse_cell_type)) / 2.5, "cm"),
    height = grid::unit(length(gene) / 2.5, "cm")
  ),
  width = 20, height = 16, units = "cm"
)
