import scanpy as sc
import os

os.makedirs('results', exist_ok=True)

def load_all_samples():
    samples = {
        '40nm': 'data/filtered_feature_bc_matrix.h5ad',
        '200nm': 'data/filtered_feature_bc_matrix_Sample2.h5ad',
        '40+200nm': 'data/filtered_feature_bc_matrix_Sample3.h5ad',
        'control': 'data/filtered_feature_bc_matrix_Sample4.h5ad'
    }
    adatas = []
    for name, path in samples.items():
        adata = sc.read_h5ad(path)
        adata.obs_names = [f"{name}_{i}" for i in adata.obs_names]
        adata.obs['sample_id'] = name
        adata.obs['condition'] = 'exposed' if name != 'control' else 'control'
        adata.var_names_make_unique()
        adatas.append(adata)
    combined_adata = sc.concat(adatas, label='sample_batch', join='outer', index_unique='_')
    return combined_adata

