import scanpy as sc
import gc
import matplotlib.pyplot as plt

def filter_data(adata, min_genes, max_genes, max_mt_pct):
    print(f'Stats BEFORE filtering: {adata.n_obs} cells, {adata.n_vars} genes')
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)
    adata = adata[adata.obs.pct_counts_mt < max_mt_pct, :].copy()
    print(f"Stats AFTER filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata

def normalize_and_select_genes(adata, n_top_genes=2500):
    print("Step 1: Saving raw counts...")
    # .raw is a special slot in AnnData. It preserves the original matrix 
    # while we perform math on the main matrix (adata.X).
    adata.raw = adata 

    print("Step 2: Normalizing total counts...")
    # Standardize library size to 10,000 counts per cell.
    # This removes the bias where some cells just have 'more' data than others.
    sc.pp.normalize_total(adata, target_sum=1e4)

    print("Step 3: Log-transforming...")
    # Apply log(x+1). This squashes the scale so that highly expressed genes 
    # don't overwhelm the statistical model.
    sc.pp.log1p(adata)

    adata.layers['norm_log'] = adata.X.copy()

    print("Step 4: Identifying Highly Variable Genes (HVGs)...")
    # We identify genes that vary significantly between our 4 samples.
    # We typically pick the top 2000 genes to focus on the 'signal' and ignore 'noise'.
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat', batch_key='sample_id')

    mt_genes = adata.var_names.str.startswith('MT-')
    ribo_genes = adata.var_names.str.startswith(('RPS', 'RPL'))

    adata.var['highly_variable'] = adata.var['highly_variable'] & ~mt_genes & ~ribo_genes
    print(f"Number of HVGs after excluding MT/Ribo: {adata.var.highly_variable.sum()}")

    adata_subset = adata[:, adata.var.highly_variable].copy()

    adata_subset.write('data/task2_ready.h5ad')

    gc.collect()

    # Visualizing HVGs for your PowerPoint
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig('results/highly_variable_genes.png')
    plt.close()

    return adata_subset