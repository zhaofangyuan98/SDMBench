import scanpy as sc
import SpaGCN as spg
from sklearn import metrics
import os
import argparse
parser = argparse.ArgumentParser()
    
# =========================== args ===============================
parser.add_argument( '--path', type=str)
args = parser.parse_args() 



adata = sc.read_h5ad(f'{args.path}')
print(adata)

adata.obs['array_row']=adata.obsm['spatial'][:,0]
adata.obs['array_col']=adata.obsm['spatial'][:,1]

x_array=adata.obs["array_row"].tolist()
y_array=adata.obs["array_col"].tolist()
adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)

for i in range(1,11):

    #Do cluster refinement
    #shape="hexagon" for Visium data, "square" for ST data.
    adj_2d_copy=adj_2d.copy()
    adata_copy=adata.copy()
    refined_pred=spg.refine(sample_id=adata_copy.obs.index.tolist(), pred=adata_copy.obs[f"pred_{i}"].tolist(), dis=adj_2d_copy, shape="hexagon")
    
    adata.obs[f"pred_refined_{i}"]=refined_pred
    adata.obs[f"pred_refined_{i}"]=adata.obs[f"pred_refined_{i}"].astype('category')
print('after:')
print(adata)
adata.write_h5ad(f'{args.path}')
