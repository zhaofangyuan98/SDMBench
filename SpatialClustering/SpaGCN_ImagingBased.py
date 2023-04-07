import scanpy as sc
from sklearn import metrics

import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import anndata as ad
import cv2


datadir = '../Data'#dataset path
memory=[1,2,3,4,5,6,7,8,9,10]
during_time=[1,2,3,4,5,6,7,8,9,10]
ari_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/osmfish.h5ad')

for j in range(10):

    import tracemalloc
    import time
 
    tracemalloc.start()  
    start_time=time.time()

    adata = sc.read_h5ad(f'{datadir}/osmfish.h5ad')
    
    adata.var_names=adata.var.index.astype("str")
    adata.var_names_make_unique()
    spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    #Normalize and take log for UMI
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    adata.obs['x_pixel']=adata.obsm['spatial'][:,0]
    adata.obs['y_pixel']=adata.obsm['spatial'][:,1]

    x_pixel=adata.obs["x_pixel"].tolist()
    y_pixel=adata.obs["y_pixel"].tolist()

    #Calculate adjacent matrix
    s=1
    b=49
    
    #If histlogy image is not available, SpaGCN can calculate the adjacent matrix using the fnction below
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
    np.savetxt(f'./adj.csv', adj, delimiter=',')
    #spatial domain detection using SpaGCN

    #expression data preprocessing
    adj=np.loadtxt(f'./adj.csv', delimiter=',')
    
    #set hyper-parameters
    p=0.5 
    #Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    n_clusters=11
    #Set seed
    r_seed=t_seed=n_seed=100
    #Seaech for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

    #run SpaGCN
    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()

    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')

    from sklearn import metrics
    obs_df = adata.obs.dropna()
    ari = metrics.adjusted_rand_score(obs_df['pred'], obs_df['Region'])
   
    ari_list[j]=ari

    end_time=time.time()
    during=end_time-start_time
    size, peak = tracemalloc.get_traced_memory()

    tracemalloc.stop()

    memory[j]=peak /1024/1024
    during_time[j]=during
    

    print('ARI score: {:.3f}'.format(ari_list[j]))
    print('memory blocks peak:{:>10.4f} MB'.format(memory[j]))
    print('time: {:.4f} s'.format(during_time[j]))


    new_data.obs['pred_{}'.format(j+1)]=adata.obs['pred']

result_path='./'
new_data.uns['time']=during_time
new_data.uns['memory']=memory
new_data.uns['ari']=ari_list
new_data.write(f'{result_path}/SpaGCN_osmfish.h5ad')


