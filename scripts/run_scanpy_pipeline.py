#!/usr/bin/env python3
"""Minimal Scanpy pipeline for single-cell neuroinformatics demos.

Accepts either:
1. .h5ad AnnData file
2. 10x mtx directory readable by scanpy.read_10x_mtx

Outputs QC plots, UMAP/clustering plots, marker tables, a processed h5ad,
and a Markdown summary suitable for GitHub.
"""

import argparse
from pathlib import Path
import scanpy as sc
import pandas as pd


def read_input(path: str):
    p = Path(path)
    if p.is_file() and p.suffix == ".h5ad":
        return sc.read_h5ad(p)
    if p.is_dir():
        return sc.read_10x_mtx(p, var_names="gene_symbols", cache=False)
    raise ValueError(f"Unsupported input: {path}. Use .h5ad or 10x mtx directory.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out_h5ad", required=True)
    parser.add_argument("--figures", required=True)
    parser.add_argument("--tables", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--min_genes", type=int, default=200)
    parser.add_argument("--max_pct_mt", type=float, default=20)
    parser.add_argument("--n_top_genes", type=int, default=2000)
    parser.add_argument("--resolution", type=float, default=0.5)
    args = parser.parse_args()

    fig_dir = Path(args.figures)
    table_dir = Path(args.tables)
    report_path = Path(args.report)
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = str(fig_dir)
    sc.settings.autoshow = False

    adata = read_input(args.input)
    adata.var_names_make_unique()

    # Basic mitochondrial QC. Human genes usually start with MT-, mouse with mt-.
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    n_cells_initial = adata.n_obs
    n_genes_initial = adata.n_vars

    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True, save="_qc.png")

    adata = adata[adata.obs["n_genes_by_counts"] >= args.min_genes].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= args.max_pct_mt].copy()

    # Standard lightweight Scanpy workflow.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.resolution, key_added="cluster")

    sc.pl.umap(adata, color=["cluster", "pct_counts_mt", "total_counts"], save="_clusters_qc.png")

    # Marker genes per cluster.
    sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(table_dir / "cluster_marker_genes.csv", index=False)

    # Optional canonical brain marker quick-look. Missing genes are skipped.
    marker_genes = ["RBFOX3", "SNAP25", "GFAP", "AQP4", "MBP", "PLP1", "CX3CR1", "P2RY12", "PDGFRA", "VCAN"]
    present = [g for g in marker_genes if g in adata.var_names]
    if present:
        sc.pl.dotplot(adata, present, groupby="cluster", save="_brain_markers.png")

    adata.write(args.out_h5ad)

    summary = f"""# Pipeline summary\n\n## Input\n- Cells before filtering: {n_cells_initial}\n- Genes before filtering: {n_genes_initial}\n\n## QC thresholds\n- Minimum genes per cell: {args.min_genes}\n- Maximum mitochondrial percent: {args.max_pct_mt}\n\n## Output\n- Cells after filtering: {adata.n_obs}\n- Highly variable genes retained: {adata.n_vars}\n- Clusters detected: {adata.obs['cluster'].nunique()}\n\n## Files generated\n- `ad_processed.h5ad`\n- `figures/` QC and UMAP plots\n- `tables/cluster_marker_genes.csv`\n\n## Biological interpretation\nUse the marker gene table and canonical brain markers to annotate clusters as neuronal, astrocytic, oligodendrocyte, microglial, OPC, or vascular-associated populations where supported by expression evidence.\n"""
    report_path.write_text(summary)


if __name__ == "__main__":
    main()
