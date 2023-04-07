import scanpy as sc
from sklearn import metrics

datadir = '../Data'#dataset path
memory=[1,2,3,4,5,6,7,8,9,10]
during_time=[1,2,3,4,5,6,7,8,9,10]
ari_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/151673.h5ad')


from sklearn import mixture

# def mclust_P(adata, num_cluster, used_obsm='STAGATE', modelNames='EEE'):
    
#     from sklearn import mixture
#     np.random.seed(2020)
#     g = mixture.GaussianMixture(n_components=num_cluster, covariance_type='diag')
#     res = g.fit_predict(adata.obsm[used_obsm])
#     adata.obs['mclust'] = res
#     return adata

for i in range(10):

    import tracemalloc
    import time
 
    tracemalloc.start()
    start_time=time.time()

    
    adata = sc.read_h5ad(f'{datadir}/151673.h5ad')
    print(adata)
    
    
    import warnings
    warnings.filterwarnings("ignore")

    import pandas as pd
    import numpy as np
    import scanpy as sc
    import matplotlib.pyplot as plt
    import os
    import sys

    from sklearn.metrics.cluster import adjusted_rand_score

    import STAGATE_pyG

    #Normalization
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    #Constructing the spatial network
    STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
    STAGATE_pyG.Stats_Spatial_Net(adata)

    #Running STAGATE
    adata = STAGATE_pyG.train_STAGATE(adata)

    sc.pp.neighbors(adata, use_rep='STAGATE')
    sc.tl.umap(adata)
    adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=7)
    # adata = mclust_P(adata, used_obsm='STAGATE', num_cluster=7))
    
    obs_df = adata.obs.dropna()
    ari = adjusted_rand_score(obs_df['mclust'], obs_df['Region'])
    ari_list[i]=ari

    end_time=time.time()
    during=end_time-start_time

    size, peak = tracemalloc.get_traced_memory()

    tracemalloc.stop()
    
    memory[i]=peak /1024/1024
    during_time[i]=during
    
    print('memory blocks peak:{:>10.4f} MB'.format(memory[i]))
    print('time: {:.4f} s'.format(during_time[i]))
    print('ARI:{}'.format(ari_list[i]))

    new_data.obs['pred_{}'.format(i+1)]=adata.obs['mclust']

result_path='./'
new_data.uns['time']=during_time
new_data.uns['memory']=memory
new_data.uns['ari']=ari_list
new_data.write(f'{result_path}/STAGATE_151673.h5ad')


