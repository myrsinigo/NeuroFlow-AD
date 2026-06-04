import pandas as pd
import scanpy as sc

markers = pd.read_csv(
    "results_sea_ad/tables/cluster_marker_genes.csv"
)

adata = sc.read_h5ad(
    "results_sea_ad/ad_processed.h5ad"
)

gene_map = adata.var["feature_name"].to_dict()

markers["gene_symbol"] = markers["names"].map(gene_map)

markers.to_csv(
    "results_sea_ad/tables/cluster_marker_genes_annotated.csv",
    index=False
)

print(markers.head(20))
