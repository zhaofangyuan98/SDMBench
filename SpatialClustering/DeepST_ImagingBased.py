import scanpy as sc
import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp

from DeepST import DeepST

datadir = '../Data'#dataset path
memory=[1,2,3,4,5,6,7,8,9,10]
during_time=[1,2,3,4,5,6,7,8,9,10]
ari_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/osmfish.h5ad')

#Run 10 repeated experiments
for i in range(10):

    import tracemalloc
    import time
 
    tracemalloc.start()
    start_time=time.time()

    
    adata = sc.read_h5ad(f'{datadir}/osmfish.h5ad')
    print(adata)
    

    # Run device
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    # the location of R, which is necessary for mclust algorithm. Please replace it with local R installation path
    os.environ['R_HOME'] ='/home/workspace/zhaofangyuan/anaconda3/envs/DeepST/lib/R'

    # the number of clusters
    n_clusters = 11

    adata.var_names_make_unique()

    # define and train model
    model = DeepST.DeepST(adata, device=device)
    adata = model.train_DeepST()

    # set radius to specify the number of neighbors considered during refinement
    radius = 50

    # clustering
    from DeepST.utils import clustering
    clustering(adata, n_clusters, radius=radius, refinement=False) #For DLPFC dataset, we use optional refinement step.

    obs_df = adata.obs.dropna()
    # calculate metric ARI
    ari = metrics.adjusted_rand_score(obs_df['domain'], obs_df['Region'])

    print('ARI:', ari)


    ari_list[i]=ari
    print('ARI:{}'.format(ari_list[i]))

    end_time=time.time()
    during=end_time-start_time
  
    size, peak = tracemalloc.get_traced_memory()

    tracemalloc.stop()
    
    memory[i]=peak /1024/1024
    during_time[i]=during
    

    print('memory blocks peak:{:>10.4f} MB'.format(memory[i]))
    print('time: {:.4f} s'.format(during_time[i]))
    
    domain=adata.obs['domain'].tolist()
    new_data.obs['pred_{}'.format(i+1)]=domain
    new_data.obs['pred_{}'.format(i+1)]=new_data.obs['pred_{}'.format(i+1)].astype('category')
    


result_path='./'
new_data.uns['time']=during_time
new_data.uns['memory']=memory
new_data.uns['ari']=ari_list
new_data.write(f'{result_path}/DeepST_osmfish.h5ad')



