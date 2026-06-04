import scanpy as sc
import numpy as np

np.random.seed(42)

input_file = "data/sea_ad.h5ad"
output_file = "data/sea_ad_subset_5000.h5ad"

adata = sc.read_h5ad(input_file, backed="r")

print(adata)
print("Available disease labels:")
print(adata.obs["disease"].value_counts())

n = min(5000, adata.n_obs)
idx = np.random.choice(adata.n_obs, size=n, replace=False)

subset = adata[idx, :].to_memory()
subset.write_h5ad(output_file)

print(f"Saved {output_file}")
print(subset)
