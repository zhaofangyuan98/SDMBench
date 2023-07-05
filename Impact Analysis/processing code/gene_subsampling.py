'''
Note that each experiment requires a re-run of the gene sub-sampling.
'''


import anndata
import numpy as np
import math

# coding:utf-8
import argparse
parser = argparse.ArgumentParser()
    
# =========================== args ===============================
parser.add_argument( '--name', type=str)
parser.add_argument( '--datadir', type=str)
args = parser.parse_args()

# read h5ad file
adata = anndata.read_h5ad(f'{args.datadir}/{args.name}.h5ad')

#gene sub-sampling
for i in range(5):
    rate=round(0.2*i+0.1,1)#rate is 10%，30%，50%，70%，90%
    print(rate)

    n_genes = math.ceil(rate*adata.n_vars)  # sub-sampling gene number
    gene_indices = np.random.choice(range(adata.n_vars), size=n_genes, replace=False)

    adata_downsampled = adata[:, gene_indices]

    print(adata_downsampled)
    adata_downsampled.write_h5ad(f'/home/workspace2/zhaofangyuan/data_h5ad/with_groundtruth/gene_downsample/{args.name}_{rate}.h5ad')