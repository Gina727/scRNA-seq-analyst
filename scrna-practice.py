import scanpy as sc
#pip install scvi-tools
#!pip install leidenalg
import scvi
import pandas as pd
import numpy as np
import seaborn as sns


adata = sc.read_csv('Desktop/GSM4285803_scRNA_RawCounts.csv')

def preprocessing(csv_path):
    adata = sc.read_csv(csv_path).T
    sc.pp.filter_genes(adata, min_cells = 10)
    sc.pp.highly_variable_genes(adata, n_top_genes = 2000, subset = True, flavor = 'seurat_v3') #select highly variable genes
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    df = solo.predict()
    df['prediction'] = solo.predict(soft = False)
    df.index = df.index.map(lambda x: x[:-2])
    df['dif'] = df.doublet - df.singlet
    doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]
    
    adata = sc.read_csv(csv_path).T
    adata.obs['Sample'] = csv_path.split('_')[2] #'raw_counts/GSM5226574_C51ctr_raw_counts.csv'
    
    adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]
    
    
    sc.pp.filter_cells(adata, min_genes=200) #get rid of cells with fewer than 200 genes
    #sc.pp.filter_genes(adata, min_cells=3) #get rid of genes that are found in fewer than 3 cells
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]
    adata = adata[adata.obs.pct_counts_mt < 20]
    adata = adata[adata.obs.pct_counts_ribo < 2]

    return adata

adata = preprocessing(csv_path)

# sc.pp.highly_variable_genes(adata, n_top_genes = 2000)
# sc.pl.highly_variable_genes(adata)
# adata = adata[:, adata.var.highly_variable]
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt', 'pct_counts_ribo'])
# sc.pp.scale(adata, max_value=10)
# sc.tl.pca(adata, svd_solver='arpack')
# sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50)
# sc.pp.neighbors(adata, n_pcs = 30) #choosing pcs 30
# sc.tl.umap(adata)
# sc.pl.umap(adata)
# sc.tl.leiden(adata, resolution = 0.5)
# #adata.obs
# sc.pl.umap(adata, color=['leiden'])

import os
sc.pp.filter_genes(adata, min_cells = 10)
from scipy.sparse import csr_matrix
adata.X = csr_matrix(adata.X)
adata.write_h5ad('combined.h5ad')

# adata.obs.groupby('Sample').count() #if multiple samples
# sc.pp.filter_genes(adata, min_cells = 100)
# adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
adata.raw = adata
scvi.model.SCVI.setup_anndata(adata, layer = "counts",
                             categorical_covariate_keys=["Sample"],
                             continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'])
model = scvi.model.SCVI(adata)
model.train() #may take a while without GPU
adata.obsm['X_scVI'] = model.get_latent_representation()
adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)
sc.pp.neighbors(adata, use_rep = 'X_scVI')
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.5)

sc.pl.umap(adata, color = ['leiden'], frameon = False)
# sc.pl.umap(adata, color = ['leiden', 'Sample'], frameon = False) #if multiple samples
adata.write_h5ad('integrated.h5ad')

# Average expression file (of each cluster)

# Cell annotation (matching the clusters with cell types)
sc.tl.leiden(adata, resolution = 1)
sc.tl.rank_genes_groups(adata, 'leiden')
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]
markers_scvi = model.differential_expression(groupby = 'leiden')
markers_scvi = markers_scvi[(markers_scvi['is_de_fdr_0.05']) & (markers_scvi.lfc_mean > .5)]
sc.pl.umap(adata, color = ['leiden'], frameon = False, legend_loc = "on data")
#...

# Differential expression graph





