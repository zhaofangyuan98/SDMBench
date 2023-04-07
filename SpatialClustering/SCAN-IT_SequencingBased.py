import scanpy as sc
from sklearn import metrics

#pre-processing
import gc
import scanit
import torch
import random
import scanpy as sc
import pandas as pd
import anndata
import numpy as np
from scipy import sparse
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
from sklearn.cluster import SpectralClustering, KMeans
import matplotlib.pyplot as plt
import stlearn as st
from pathlib import Path

datadir = '../Data'#dataset path
memory=[1,2,3,4,5,6,7,8,9,10]
during_time=[1,2,3,4,5,6,7,8,9,10]
ari_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/151673.h5ad')


def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    '''
        arg1(adata)[AnnData matrix]
        arg2(fixed_clus_count)[int]
        
        return:
            resolution[int]
    '''
    for res in sorted(list(np.arange(0.15, 5, increment)), reverse=True):
    # for res in sorted(list(np.arange(0.2, 2.5, increment)), reverse=True):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        # print(res)
        # print(count_unique_leiden)
        if count_unique_leiden == fixed_clus_count:
            break
    return res

for i in range(10):
    import tracemalloc
    import time
 
    tracemalloc.start()
    start_time=time.time()

    adata = sc.read_h5ad(f'{datadir}/151673.h5ad')
    print(adata)

    #SOMDE
    from somde import SomNode
    pts = adata.obsm['spatial']
    df_sp = pd.DataFrame(data=adata.X, columns=list(adata.var_names))
    
    som = SomNode(pts, 5)

    ndf,ninfo = som.mtx(df_sp.T)
    nres = som.norm()
    result, SVnum =som.run()
    result.to_csv(f'./somde_result.csv')

    adata = sc.read_h5ad(f'{datadir}/151673.h5ad')
    n_sv_genes = 3000
    adata_sp = adata.copy()
    sc.pp.normalize_total(adata_sp)
    df_somde = pd.read_csv(f'./somde_result.csv')
    sv_genes = list( df_somde['g'].values[:n_sv_genes] )
    adata_sp = adata_sp[:, sv_genes]
    sc.pp.log1p(adata_sp)
    sc.pp.scale(adata_sp)

    scanit.tl.spatial_graph(adata_sp, method='alpha shape', alpha_n_layer=2, knn_n_neighbors=5)
    scanit.tl.spatial_representation(adata_sp, n_h=30, n_epoch=2000, lr=0.001, device='cuda', n_consensus=1, projection='mds', 
                python_seed=0, torch_seed=0, numpy_seed=0)
    
    n_clusters=7

    sc.pp.neighbors(adata_sp, use_rep='X_scanit', n_neighbors=15)
    eval_resolution = res_search_fixed_clus(adata_sp, n_clusters)
    sc.tl.leiden(adata_sp, key_added="scanit_leiden", resolution=eval_resolution)

    from sklearn import metrics
    obs_df = adata_sp.obs.dropna()
    ari = adjusted_rand_score(obs_df['scanit_leiden'], obs_df['Region'])
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

    new_data.obs['pred_{}'.format(i+1)]=adata_sp.obs['scanit_leiden']

result_path='./'
new_data.uns['time']=during_time
new_data.uns['memory']=memory
new_data.uns['ari']=ari_list
new_data.write(f'{result_path}/SCAN-IT_151673.h5ad')


