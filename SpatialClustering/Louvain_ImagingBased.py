import scanpy as sc
from sklearn import metrics

datadir = '../Data'#dataset path
memory=[1,2,3,4,5,6,7,8,9,10]
during_time=[1,2,3,4,5,6,7,8,9,10]
ari_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/osmfish.h5ad')

def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):

    for res in sorted(list(np.arange(0.2, 2.5, increment)), reverse=True):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        if count_unique_leiden == fixed_clus_count:
            break
   
    return res


#Run 10 repeated experiments
for i in range(10):

    import tracemalloc
    import time
 
    tracemalloc.start()  
    start_time=time.time()

    
    adata = sc.read_h5ad(f'{datadir}/osmfish.h5ad')

    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
   
    n_clusters=11
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50) # 20
    eval_resolution = res_search_fixed_clus(adata, n_clusters)
        
    sc.tl.louvain(adata, key_added="louvain", resolution=eval_resolution)

    obs_df = adata.obs.dropna()
    ari = metrics.adjusted_rand_score(obs_df['Region'], obs_df['louvain'])
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

    new_data.obs['pred_{}'.format(i+1)]=adata.obs['louvain']

result_path='./'
new_data.uns['time']=during_time
new_data.uns['memory']=memory
new_data.uns['ari']=ari_list
new_data.write(f'{result_path}/Louvain_osmfish.h5ad')

