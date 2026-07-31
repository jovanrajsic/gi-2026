import scanpy as sc
import matplotlib.pyplot as plt
import scanpy.external as sce

def run_pca(adata):
    # 1. Scale the data (Only on the HVGs to save RAM!)
    # max_value=10 prevents highly expressed outliers from dominating
    sc.pp.scale(adata, max_value=10)

    # 2. Run PCA
    # svd_solver='arpack' is the most memory-efficient for 8GB RAM
    sc.tl.pca(
        adata,
        n_comps=50,                # Calculate top 50 PCs
        use_highly_variable=True,   # Safety check (already sliced, but good practice)
        svd_solver="arpack",        # Fast/Stable for your 33k cells
        random_state=0              # Makes your results same every time
    )

    # 3. Visualization: The Elbow Plot (Variance Ratio)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
    plt.savefig('results/pca_elbow_plot_log.png')
    plt.show()

    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=False, show=False)
    plt.savefig('results/pca_elbow_plot.png')
    plt.show()

    sc.pl.pca(adata, color='sample_id', show=False)
    plt.savefig('results/pca_scatter.png')
    plt.show()

    return adata

def run_clustering(adata, n_pcs=30, resolution=0.4):
    # batch_key is the column in adata.obs that identifies the 4 samples
    sce.pp.bbknn(adata, batch_key='sample_id', n_pcs=n_pcs)

    sc.tl.umap(
        adata, 
        min_dist=0.1,    # Keep it at 0.4 to handle the 'dropout' noise
        spread=1.0,      # Good for readability
        random_state=0   # For reproducibility
    )

    sc.tl.leiden(
        adata,
        key_added="leiden",
        resolution=resolution, # Larger resolution - more clusters
        random_state=0
    )

    # Plot it!
    sc.pl.umap(adata, color=['sample_id', 'leiden', 'pct_counts_mt'], cmap='turbo', show=False)
    plt.savefig('results/umap.png')
    plt.show()

    return adata