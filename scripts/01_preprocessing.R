# ---------------------------
# Title: 01_preprocessing — QC, normalization, Harmony, clustering, UMAP, coarse annotation
# Purpose:
#   - Works for either tumor or lymph node (LN) Seurat objects.
#   - QC thresholds and PCA/Harmony dims are identical for tumor and LN.
#   - Only clustering resolution and cluster→coarse annotation differ by dataset.
# Inputs:
#   - Raw Seurat object in memory: `seu_raw` (loaded upstream)
#   - annotations/tumor_coarse_cluster_mapping.tsv   (cluster, coarse_cell_type)
#   - annotations/ln_coarse_cluster_mapping.tsv      (cluster, coarse_cell_type)
# Outputs:
#   - ../data/Tumor_scRNA-seq_seurat.rds   (if DATASET == "tumor")
#   - ../data/LN_scRNA-seq_seurat.rds      (if DATASET == "ln")
# Notes:
#   - Doublets flagged by Scrublet were removed upstream.
#   - Cells with high skeletal-muscle marker expression and globally high ribosome expression
#     were excluded prior to this step.
# ---------------------------

# Packages
library(Seurat)
library(dplyr)
library(harmony)

# ---------------------------
# 0) Dataset switch (only things that differ by tissue)
# ---------------------------
DATASET <- "tumor"   # set to "tumor" or "ln"

RES_CLUSTER_TUMOR <- 1.0
RES_CLUSTER_LN    <- 0.8     # choose your LN resolution here

res_cluster <- if (DATASET == "tumor") RES_CLUSTER_TUMOR else RES_CLUSTER_LN
anno_path   <- if (DATASET == "tumor") {
  "annotations/tumor_coarse_cluster_mapping.tsv"
} else {
  "annotations/ln_coarse_cluster_mapping.tsv"
}
out_rds <- if (DATASET == "tumor") {
  "../data/Tumor_scRNA-seq_seurat.rds"
} else {
  "../data/LN_scRNA-seq_seurat.rds"
}

# ---------------------------
# 1) QC metrics (shared by tumor and LN)
# ---------------------------
# Mouse mitochondrial genes often start with "mt-"; adjust if needed (e.g., "^MT-")
seu_raw[["percent.mt"]] <- PercentageFeatureSet(seu_raw, pattern = "^mt-")
seu_raw[["percent.rb"]] <- PercentageFeatureSet(seu_raw, pattern = "^Rp[sl]")

# ---------------------------
# 2) Cell filtering (shared thresholds)
# ---------------------------
seu_qc <- subset(
  seu_raw,
  subset =
    nFeature_RNA > 1000 &
    nFeature_RNA < 8000 &
    nCount_RNA   > 2000 &
    nCount_RNA   < 50000 &
    percent.mt   < 5
)

# ---------------------------
# 3) Normalize → HVG → Scale (regress mt/ribo) → PCA (shared)
# ---------------------------
seu_norm <- seu_qc %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c("percent.mt", "percent.rb")) %>%
  RunPCA(npcs = 50)

# ---------------------------
# 4) Harmony batch correction (shared)
# ---------------------------
seu_norm <- RunHarmony(
  object        = seu_norm,
  group.by.vars = "orig.ident",
  assay.use     = "RNA"
)

# ---------------------------
# 5) Graph / clustering / UMAP on Harmony dims (shared dims, dataset-specific resolution)
# ---------------------------
dims_use <- 1:50

seu_norm <- FindNeighbors(seu_norm, reduction = "harmony", dims = dims_use)
seu_norm <- FindClusters(seu_norm, resolution = res_cluster)
seu_norm <- RunUMAP(seu_norm, reduction = "harmony", dims = dims_use, min.dist = 0.3)

# ---------------------------
# 6) Coarse annotation: cluster → coarse_cell_type (dataset-specific mapping)
#    Mapping files are simple two-column TSVs:
#      cluster   coarse_cell_type
#      0         Tumor cell
#      1         T/NK cell
#      ...
# ---------------------------
cluster_map <- read.table(anno_path, header = TRUE, sep = "\t", check.names = FALSE)

meta <- seu_norm@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(cluster = as.integer(as.character(seurat_clusters))) %>%
  left_join(cluster_map, by = c("cluster" = "cluster"))

# fallback if any cluster is not mapped
if (any(is.na(meta$coarse_cell_type))) {
  meta$coarse_cell_type[is.na(meta$coarse_cell_type)] <- "Other"
}

seu_norm@meta.data <- meta %>% tibble::column_to_rownames("cell_id")

# ---------------------------
# 7) Save
# ---------------------------
saveRDS(seu_norm, file = out_rds)