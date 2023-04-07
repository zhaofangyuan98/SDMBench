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
seed_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/151673.h5ad')

for j in range(10):

    import tracemalloc
    import time
 
    tracemalloc.start() 
    start_time=time.time()

    adata = sc.read_h5ad(f'{datadir}/151673.h5ad')
    print(adata)
    
    
    spatial=pd.read_csv(f"./positions.txt",sep=",",header=None,na_filter=False,index_col=0) 
    adata.obs["x1"]=spatial[1]
    adata.obs["x2"]=spatial[2]
    adata.obs["x3"]=spatial[3]
    adata.obs["x4"]=spatial[4]
    adata.obs["x5"]=spatial[5]
    adata.obs["x_array"]=adata.obs["x2"]
    adata.obs["y_array"]=adata.obs["x3"]
    adata.obs["x_pixel"]=adata.obs["x4"]
    adata.obs["y_pixel"]=adata.obs["x5"]
    obsm_spatial=spatial.loc[adata.obs_names].loc[:,[4,5]]
    adata.obsm['spatial']=obsm_spatial.values
    adata.obsm['spatial']
    
    #Select captured samples
    adata=adata[adata.obs["x1"]==1]
    adata.var_names=[i.upper() for i in list(adata.var_names)]#将字符串的小写字母变成大写字母
    adata.var["genename"]=adata.var.index.astype("str")
    
    #expression data preprocessing
    
    adata.var_names=adata.var.index.astype("str")
    adata.var_names_make_unique()
    spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    #Normalize and take log for UMI
    sc.pp.normalize_per_cell(adata)
    print(adata)
    sc.pp.log1p(adata)

    #Read in hitology image
    img=cv2.imread(f"./image.tif")

    #integrate gene expression and histilogy into a Graph

    #Set coordinates
    x_array=adata.obs["x_array"].tolist()
    y_array=adata.obs["y_array"].tolist()
    x_pixel=adata.obs["x_pixel"].tolist()
    y_pixel=adata.obs["y_pixel"].tolist()

    #Test coordinates on the image
    img_new=img.copy()
    for i in range(len(x_pixel)):
        x=x_pixel[i]
        y=y_pixel[i]
        img_new[int(x-20):int(x+20), int(y-20):int(y+20),:]=0
    cv2.imwrite(f'./151673_map.jpg', img_new)

    #Calculate adjacent matrix
    s=1
    b=49
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
    np.savetxt(f'./adj.csv', adj, delimiter=',')#3639*3639的邻接矩阵

    #spatial domain detection using SpaGCN

    adj=np.loadtxt(f'./adj.csv', delimiter=',')

    #set hyper-parameters
    p=0.5 
    #Find the l value given p
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    n_clusters=7
    #Set seed
    r_seed=t_seed=n_seed=seed_list[j]
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
    obs_df=adata.obs.dropna()
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
new_data.write(f'{result_path}/SpaGCN(HE)_151673.h5ad')


