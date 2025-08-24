#!/usr/bin/env python3
# ---------------------------
# Title: 01_preprocessing_scrublet.py
# Purpose:
#   - Detect doublets in single-cell RNA-seq data using Scrublet
#   - Save histogram and UMAP plots
#   - Export per-cell doublet scores and predictions
# Inputs:
#   - 10x Genomics matrices: matrix.mtx(.gz), barcodes.tsv(.gz)
# Outputs:
#   - <sample>.histogram.pdf
#   - <sample>.umap.pdf
#   - <sample>.doublet.tsv
# Notes:
#   - Reproducibility ensured with fixed SEED
# ---------------------------

import gzip
import numpy as np
import pandas as pd
import scipy.io
import scrublet as scr
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------
# Parameters
# ---------------------------
INPUT_DIR   = Path("/path/to/10x_root")    # root directory containing sample folders
SAMPLES     = ["Sample1", "Sample2"]       # edit with your sample names
OUTPUT_DIR  = Path("./scrublet_output")

EXPECTED_RATE_PER_1K = 0.008   # expected doublet rate per 1000 cells
MIN_COUNTS           = 3
MIN_CELLS            = 3
MIN_GENE_PCTL        = 85
N_PCS                = 30

UMAP_N_NEIGHBORS     = 10
UMAP_MIN_DIST        = 0.3

# ---------------------------
# Helpers
# ---------------------------
def read_barcodes(path: Path) -> pd.DataFrame:
    """Read 10x barcodes.tsv(.gz) and return DataFrame with column 'barcode'."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            df = pd.read_csv(fh, header=None, names=["barcode"])
    else:
        df = pd.read_csv(path, header=None, names=["barcode"])
    return df

def read_matrix(path: Path):
    """Read 10x matrix.mtx(.gz) and return cell × gene sparse matrix (CSC)."""
    mtx = scipy.io.mmread(str(path))
    return mtx.T.tocsc()

def run_scrublet(sample: str):
    """Run Scrublet for one sample and export results."""
    sample_dir = INPUT_DIR / sample

    # locate files
    mtx = sample_dir / "matrix.mtx"
    mtx_gz = sample_dir / "matrix.mtx.gz"
    bc = sample_dir / "barcodes.tsv"
    bc_gz = sample_dir / "barcodes.tsv.gz"

    mtx_path = mtx_gz if mtx_gz.exists() else mtx
    bc_path  = bc_gz if bc_gz.exists() else bc
    if not mtx_path.exists() or not bc_path.exists():
        raise FileNotFoundError(f"Missing 10x files for sample {sample}")

    # load
    counts   = read_matrix(mtx_path)
    barcodes = read_barcodes(bc_path)

    # expected doublet rate
    expected_rate = EXPECTED_RATE_PER_1K * (counts.shape[0] / 1000.0)
    scrub = scr.Scrublet(counts_matrix=counts, expected_doublet_rate=expected_rate)

    # run
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=MIN_COUNTS,
        min_cells=MIN_CELLS,
        min_gene_variability_pctl=MIN_GENE_PCTL,
        n_prin_comps=N_PCS
    )

    # histogram
    scrub.plot_histogram()
    plt.savefig(OUTPUT_DIR / f"{sample}.histogram.pdf", bbox_inches="tight")

    # UMAP embedding
    scrub.set_embedding(
        "UMAP",
        scr.get_umap(scrub.manifold_obs_, n_neighbors=UMAP_N_NEIGHBORS, min_dist=UMAP_MIN_DIST)
    )
    scrub.plot_embedding('UMAP', order_points=True)
    plt.savefig(OUTPUT_DIR / f"{sample}.UMAP.pdf", bbox_inches="tight")

    # export results
    out = barcodes.copy()
    out["doublet_score"]    = doublet_scores
    out["predicted_doublet"] = predicted_doublets.astype(bool)
    out.to_csv(OUTPUT_DIR / f"{sample}.doublet.tsv", sep="\t", index=False)

# ---------------------------
# Main
# ---------------------------
if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for s in SAMPLES:
        print(f"[Scrublet] Processing {s} ...")
        run_scrublet(s)
    print("[Scrublet] Done.")