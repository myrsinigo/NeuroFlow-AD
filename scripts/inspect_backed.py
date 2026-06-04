import scanpy as sc

adata = sc.read_h5ad("data/sea_ad.h5ad", backed="r")

print(adata)

print("\nobs columns:")
for col in list(adata.obs.columns)[:50]:
    print("-", col)

print("\nvar columns:")
for col in list(adata.var.columns)[:20]:
    print("-", col)
