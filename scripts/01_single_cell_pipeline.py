import scanpy as sc

# =========================
# 1. Load dataset
# =========================
adata = sc.datasets.pbmc3k()

# =========================
# 2. QC metrics
# =========================
adata.var['mt'] = adata.var_names.str.startswith('MT_')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# =========================
# 3. Filtering
# =========================
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.pct_counts_mt < 10, :]
sc.pp.filter_genes(adata, min_cells=3)

# =========================
# 4. Normalization
# =========================
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# =========================
# 5. HVGs
# =========================
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
adata = adata[:, adata.var.highly_variable].copy()

# =========================
# 6. PCA
# =========================
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# =========================
# 7. Neighbors + UMAP
# =========================
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# =========================
# 8. Leiden clustering
# =========================
sc.tl.leiden(adata, resolution=0.5)

# =========================
# 9. Save outputs
# =========================
adata.obs['leiden'].value_counts().to_csv('cluster_distribution.csv')
adata.write('processed_pbmc.h5ad')

print('Pipeline completed successfully')
