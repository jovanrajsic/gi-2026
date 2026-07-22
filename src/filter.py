import scanpy as sc

def filter_data(adata, min_genes, max_genes, max_mt_pct):
    print(f'Stats BEFORE filtering: {adata.n_obs} cells, {adata.n_vars} genes')
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    adata = adata[adata.obs.pct_counts_mt < max_mt_pct, :].copy()
    print(f"Stats AFTER filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata