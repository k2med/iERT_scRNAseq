# ---------------------------
# Title: Fig.S6 — CellChat Analysis and Visualization
# Purpose:
#   - Perform CellChat analysis for different treatment groups (Ctrl, LDRT, SBRT, ERT)
#   - Generate bubble plots for cell-cell interactions (Fig.S6E, Fig.S6F)
# Inputs:
#   - ../data/Tumor_scRNA-seq_seurat.rds
#   - CellChatDB.mouse (CellChat database)
# Outputs:
#   - Fig.S6E.pdf                        (Bubble plot for NK - cDC interactions)
#   - Fig.S6F.pdf                        (Bubble plot for cDC - CD8 T interactions)
#   - FigS6E_dotplot_source_data.txt      (Source data for Fig.S6E)
#   - FigS6F_dotplot_source_data.txt      (Source data for Fig.S6F)
# ---------------------------

# Packages
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CellChat)

# 0) Load data
seu_all <- readRDS("../data/Tumor_scRNA-seq_seurat.rds")

# 1) Fill NA values in 'fine_cell_type' with 'coarse_cell_type'
seu_all@meta.data[is.na(seu_all$fine_cell_type), "fine_cell_type"] <- as.vector(seu_all@meta.data[is.na(seu_all$fine_cell_type), "coarse_cell_type"])

# 2) Set active identity to 'fine_cell_type'
seu_all@active.ident <- factor(seu_all$fine_cell_type)

# 3) Perform CellChat analysis for each group
groups <- unique(seu_all$group)
for (current_group in groups) {
  seurat_subset <- subset(seu_all, group == current_group)
  
  # Create CellChat object
  cellchat <- createCellChat(object = seurat_subset, group.by = "fine_cell_type", assay = "RNA")
  
  # Load database and perform analysis
  cellchat@DB <- CellChatDB.mouse
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- smoothData(cellchat, adj = PPI.human)
  cellchat <- computeCommunProb(cellchat, population.size = TRUE, raw.use = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Save the results for each group
  saveRDS(cellchat, file = paste0(current_group, "_cellchat.rds"))
}

# 4) Load all CellChat objects for each treatment group
treatment_groups <- c("Ctrl", "LDRT", "SBRT", "ERT")
cellchat_objects <- lapply(treatment_groups, function(group) readRDS(paste0(group, "_cellchat.rds")))

names(cellchat_objects) <- treatment_groups

# 5) Define the cell types for consistent ordering
cell_types <- c('Tumor cell', 'Macrophage', 'CD4+ T cell', 'Fibroblast', 'Endothelial cell',
                "Naive/memory CD8+ T", "Precursor of exhausted CD8+ T (Tpex)", "Effector-like CD8+ T (Tef)",
                "Exhausted CD8+ T (Tex)", "ISG+ CD8+ T", "Proliferating CD8+ T", 'pDC', 'Neutrophil', 'B cell',
                'Xcl1+ NK', 'Ccl5+ NK', 'Ccr2+ NK', 'Gzmg+ NK', 'Mki67+ NK', 'Xcr1+ cDC1', 'Ccr7+ cDC1', 'Cd209a+ cDC2', 'Mki67+ cDC')

# 6) Set consistent factor levels for all CellChat objects
cellchat_objects <- lapply(cellchat_objects, function(cellchat) {
  cellchat@idents <- factor(cellchat@idents, levels = cell_types)
  return(cellchat)
})

# 7) Merge all CellChat objects for comparison
cellchat <- mergeCellChat(cellchat_objects, add.names = names(cellchat_objects))

# ---------------------------
# 8) Visualize interactions: NK vs cDC (Fig.S6E)
# ---------------------------
# Define pair interactions to visualize
pairLR.use <- data.frame(interaction_name = rev(c('XCL1_XCR1', 'CCL5_CCR5', 'CCL4_CCR5', 'CCL3_CCR5', 'CCL3_CCR1',
                                                  'CCL8_CCR5', 'CCL8_CCR2', 'CCL8_CCR1', 'CCL2_CCR2', 'CXCL10_CXCR3')))

# Create bubble plot for NK vs cDC interactions
plot <- netVisual_bubble(cellchat, sources.use = c("Xcl1+ NK", "Ccl5+ NK", "Ccr2+ NK", "Gzmg+ NK"),
                         targets.use = c("Xcr1+ cDC1", "Ccr7+ cDC1", "Cd209a+ cDC2"),
                         pairLR.use = pairLR.use, comparison = c(1, 2, 3, 4), angle.x = 45)

# Save plot as PDF
ggsave("Fig.S6E.pdf", width = 30, height = 8, units = "cm",
       egg::set_panel_size(plot, width = unit(48 * 0.4, "cm"), height = unit(10 * 0.4, "cm")))

# Save source data for Fig.S6E
dot_df <- plot$data[, c("interaction_name_2", "source", "target", "dataset", "prob", "pval")]
dot_df <- dot_df[order(dot_df$interaction_name_2, dot_df$source, dot_df$target, decreasing = c(TRUE, FALSE, FALSE)),]
dot_df$pval <- case_when((dot_df$pval == 3) ~ "p < 0.01",
                         (dot_df$pval == 2) ~ "0.01 < p < 0.05",
                         (dot_df$pval == 1) ~ "p > 0.05")
colnames(dot_df) <- c("interaction", "sender", "receiver", "group", "prob", "pval")

write.table(dot_df, "FigS6E_dotplot_source_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# 9) Visualize interactions: cDC vs CD8 T (Fig.S6F)
# ---------------------------
# Define pair interactions to visualize
pairLR.use <- data.frame(interaction_name = rev(c('H2-K1_CD8A', 'H2-K1_CD8B1', 'H2-D1_CD8A', 'H2-D1_CD8B1')))

# Create bubble plot for cDC vs CD8 T interactions
plot <- netVisual_bubble(cellchat, sources.use = c("Xcr1+ cDC1", "Ccr7+ cDC", "Cd209a+ cDC2"),
                         targets.use = c("Naive/memory CD8+ T", "Precursor of exhausted CD8+ T (Tpex)", "Effector-like CD8+ T (Tef)",
                                         "Exhausted CD8+ T (Tex)", "ISG+ CD8+ T"), pairLR.use = pairLR.use,
                         comparison = c(1, 2, 3, 4), angle.x = 45)

# Save plot as PDF
ggsave("Fig.S6F.pdf", width = 40, height = 8, units = "cm",
       egg::set_panel_size(plot, width = unit(60 * 0.4, "cm"), height = unit(4 * 0.4, "cm")))

# Save source data for Fig.S6F
dot_df <- plot$data[, c("interaction_name_2", "source", "target", "dataset", "prob", "pval")]
dot_df <- dot_df[order(dot_df$interaction_name_2, dot_df$source, dot_df$target, decreasing = c(TRUE, FALSE, FALSE)),]
dot_df$pval <- case_when((dot_df$pval == 3) ~ "p < 0.01",
                         (dot_df$pval == 2) ~ "0.01 < p < 0.05",
                         (dot_df$pval == 1) ~ "p > 0.05")
colnames(dot_df) <- c("interaction", "sender", "receiver", "group", "prob", "pval")

write.table(dot_df, "FigS6F_dotplot_source_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)
