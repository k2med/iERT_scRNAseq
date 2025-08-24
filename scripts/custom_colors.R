# ===========================
# Custom color palettes
# ===========================

# ---------------------------
# Overall UMAP palette
# Includes: major cell types + experimental conditions
# ---------------------------
overall_umap_colors <- c(
  # cell types
  "Tumor cell"        = "#7194b8",
  "CD8+ T cell"       = "#e7797a",
  "CD4+ T cell"       = "#f4a454",
  "Gamma delta T cell"= "#2a7fc3",
  "NK cell"           = "#a38cbd",
  "Macrophage"        = "#bb7e50",
  "Neutrophil"        = "#79b473",
  "cDC"               = "#d7c23c",
  "pDC"               = "#a5d38f",
  "B cell"            = "#e16712",
  "Fibroblast"        = "#90c4c0",
  "Endothelial cell"  = "#ab9697",
  
  # experimental conditions
  "ERT"               = "#118002",
  "SBRT"              = "#0000ff",
  "LDRT"              = "#fd8008",
  "Ctrl"              = "#d9dad8"
)

# ---------------------------
# NK subset palette
# ---------------------------
nk_umap_colors <- c(
  "Xcl1+ NK"  = "#dd7a80",
  "Ccl5+ NK"  = "#76bbd8",
  "Ccr2+ NK"  = "#bbb7cb",
  "Gzmg+ NK"  = "#aacb65",
  "Mki67+ NK" = "#eb9c7f"
)

# ---------------------------
# cDC subset palette
# ---------------------------
cdc_umap_colors <- c(
  "Xcr1+ cDC1"   = "#ec8c70",
  "Ccr7+ cDC1"   = "#dbd08c",
  "Cd209a+ cDC2" = "#a0d09d",
  "Mki67+ cDC"   = "#4690a9"
)

# ---------------------------
# CD8+ T subset palette
# ---------------------------
cd8t_umap_colors <- c(
  "Naive/memory CD8+ T"              = "#80b1d3",
  "Precursor of exhausted CD8+ T (Tpex)" = "#ec8c70",
  "Effector-like CD8+ T (Tef)"       = "#45998b",
  "Exhausted CD8+ T (Tex)"           = "#dfc987",
  "ISG+ CD8+ T"                      = "#d4b9da",
  "Proliferating CD8+ T"             = "#fdb462"
)

# ---------------------------
# Gene expression heatmap palette
# ---------------------------
heatmap_colors <- c(
  "black", "#006aff", "#0068fe", "#0066fc", "#0063fa", "#0061f8", "#005ff6", "#005df4",
  "#005af2", "#0058f0", "#0056ee", "#0054ec", "#0052ea", "#004fe8", "#004de6", "#004be4",
  "#0049e2", "#0046e0", "#0045de", "#0044dd", "#0042db", "#0041d9", "#0040d8", "#003fd6",
  "#003ed4", "#003dd3", "#003cd1", "#003bcf", "#003ace", "#0039cc", "#0038ca", "#0037c9",
  "#0036c7", "#0034c5", "#0033c4", "#0131ba", "#0c2eae", "#132ca1", "#172995", "#192688",
  "#1b247c", "#1c2170", "#1c1f65", "#1b1c59", "#1a1a4e", "#191743", "#171539", "#15122e",
  "#130e25", "#100a1b", "#080511", "#000000", "#121205", "#1d1e0a", "#26290e", "#30360f",
  "#3b4210", "#464f0f", "#515c0e", "#5c6a0c", "#677807", "#738601", "#7f9500", "#8ba400",
  "#97b300", "#a4c200", "#b1d200", "#bde100", "#c6ec00", "#c8ee00", "#c9ef00", "#caf000",
  "#cbf100", "#ccf300", "#cef400", "#cff500", "#d0f600", "#d1f800", "#d2f900", "#d4fa00",
  "#d5fc00", "#d6fd00", "#d7fe00", "#d8ff00", "#daff00", "#dbff00", "#ddff00", "#dfff00",
  "#e0ff00", "#e2ff00", "#e4ff00", "#e6ff00", "#e7ff00", "#e9ff00", "#ebff00", "#edff02",
  "#eeff07", "#f0ff0b", "#f2ff0e", "#f4ff12", "#f5ff14", "#f7ff17"
)
