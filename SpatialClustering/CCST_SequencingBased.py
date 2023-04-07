import scanpy as sc
from sklearn import metrics

import pandas as pd
import scanpy as sc
import numpy as np
import stlearn as st

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
from h5py import Dataset, Group

import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn import metrics
from scipy import sparse
#from sklearn.metrics import roc_curve, auc, roc_auc_score

import anndata as ad
import numpy as np
import pickle
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, ChebConv, GATConv, DeepGraphInfomax, global_mean_pool, global_max_pool  # noqa
from torch_geometric.data import Data, DataLoader
from datetime import datetime 

from sklearn.metrics.cluster import adjusted_rand_score
from CCST import get_graph, train_DGI, train_DGI, PCA_process, Kmeans_cluster

datadir = '../Data'#dataset path
memory=[1,2,3,4,5,6,7,8,9,10]
during_time=[1,2,3,4,5,6,7,8,9,10]
ari_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/151673.h5ad')


def adata_preprocess(i_adata, min_cells=3, pca_n_comps=300):
    print('===== Preprocessing Data ')
    sc.pp.filter_genes(i_adata, min_cells=min_cells)
    adata_X = sc.pp.normalize_total(i_adata, target_sum=1, exclude_highly_expressed=True, inplace=False)['X']
    adata_X = sc.pp.scale(adata_X)
    adata_X = sc.pp.pca(adata_X, n_comps=pca_n_comps)
    return adata_X



def get_adj(generated_data_fold):
    coordinates = np.load(generated_data_fold + 'coordinates.npy')
    if not os.path.exists(generated_data_fold):
        os.makedirs(generated_data_fold) 
    ############# get batch adjacent matrix
    cell_num = len(coordinates)

    ############ the distribution of distance 
    if 1:#not os.path.exists(generated_data_fold + 'distance_array.npy'):
        distance_list = []
        print ('calculating distance matrix, it takes a while')
        
        distance_list = []
        for j in range(cell_num):
            for i in range (cell_num):
                if i!=j:
                    distance_list.append(np.linalg.norm(coordinates[j]-coordinates[i]))

        distance_array = np.array(distance_list)
        #np.save(generated_data_fold + 'distance_array.npy', distance_array)
    else:
        distance_array = np.load(generated_data_fold + 'distance_array.npy')

    ###try different distance threshold, so that on average, each cell has x neighbor cells, see Tab. S1 for results
    from scipy import sparse
    import pickle
    import scipy.linalg

    for threshold in [300]:#range (210,211):#(100,400,40):
        num_big = np.where(distance_array<threshold)[0].shape[0]
        print (threshold,num_big,str(num_big/(cell_num*2))) #300 22064 2.9046866771985256
        from sklearn.metrics.pairwise import euclidean_distances

        distance_matrix = euclidean_distances(coordinates, coordinates)
        distance_matrix_threshold_I = np.zeros(distance_matrix.shape)
        distance_matrix_threshold_W = np.zeros(distance_matrix.shape)
        for i in range(distance_matrix_threshold_I.shape[0]):
            for j in range(distance_matrix_threshold_I.shape[1]):
                if distance_matrix[i,j] <= threshold and distance_matrix[i,j] > 0:
                    distance_matrix_threshold_I[i,j] = 1
                    distance_matrix_threshold_W[i,j] = distance_matrix[i,j]
            
        
        ############### get normalized sparse adjacent matrix
        distance_matrix_threshold_I_N = np.float32(distance_matrix_threshold_I) ## do not normalize adjcent matrix
        distance_matrix_threshold_I_N_crs = sparse.csr_matrix(distance_matrix_threshold_I_N)
        with open(generated_data_fold + 'Adjacent', 'wb') as fp:
            pickle.dump(distance_matrix_threshold_I_N_crs, fp)


def get_data(args):
    data_file = args.data_path + args.data_name +'/'
    with open(data_file + 'Adjacent', 'rb') as fp:
        adj_0 = pickle.load(fp)
    X_data = np.load(data_file + 'features.npy')

    num_points = X_data.shape[0]
    adj_I = np.eye(num_points)
    adj_I = sparse.csr_matrix(adj_I)
    adj = (1-args.lambda_I)*adj_0 + args.lambda_I*adj_I

    # cell_type_indeces = np.load(data_file + 'cell_types.npy')
    
    return adj_0, adj, X_data



def res_search_fixed_clus(cluster_type, adata, fixed_clus_count, increment=0.02):
    '''
        arg1(adata)[AnnData matrix]
        arg2(fixed_clus_count)[int]
        
        return:
            resolution[int]
    '''
    if cluster_type == 'leiden':
        for res in sorted(list(np.arange(0.2, 2.5, increment)), reverse=True):
            sc.tl.leiden(adata, random_state=0, resolution=res)
            count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
            if count_unique_leiden == fixed_clus_count:
                break
    elif cluster_type == 'louvain':
        for res in sorted(list(np.arange(0.2, 2.5, increment)), reverse=True):
            sc.tl.louvain(adata, random_state=0, resolution=res)
            count_unique_louvain = len(pd.DataFrame(adata.obs['louvain']).louvain.unique())
            if count_unique_louvain == fixed_clus_count:
                break
    return res

##Run 10 repeated experiments
for i in range(10):

    import tracemalloc
    import time
 
    tracemalloc.start()
    start_time=time.time()
 
    adata = sc.read_h5ad(f'{datadir}/151673.h5ad')
    print(adata)
    
    rootPath = os.path.dirname(sys.path[0])
    os.chdir(rootPath+'/CCST')
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--min_cells', type=float, default=5, help='Lowly expressed genes which appear in fewer than this number of cells will be filtered out')
    parser.add_argument( '--Dim_PCA', type=int, default=200, help='The output dimention of PCA')
    parser.add_argument( '--data_path', type=str, default='dataset/', help='The path to dataset')
    parser.add_argument( '--data_name', type=str, default='151673', help='The name of dataset')
    parser.add_argument( '--generated_data_path', type=str, default='generated_data/', help='The folder to store the generated data')

    # =========================== args ===============================
    parser.add_argument( '--data_name', type=str, default='151673', help="'MERFISH' or 'V1_Breast_Cancer_Block_A_Section_1") 
    parser.add_argument( '--lambda_I', type=float, default=0.3) #0.3 on 10x visium, 0.8 on other technology
    parser.add_argument( '--data_path', type=str, default='generated_data/', help='data path')
    parser.add_argument( '--model_path', type=str, default='model') 
    parser.add_argument( '--embedding_data_path', type=str, default='Embedding_data') 
    parser.add_argument( '--result_path', type=str, default='results') 
    parser.add_argument( '--DGI', type=int, default=1, help='run Deep Graph Infomax(DGI) model, otherwise direct load embeddings')
    parser.add_argument( '--load', type=int, default=0, help='Load pretrained DGI model')
    parser.add_argument( '--num_epoch', type=int, default=5000, help='numebr of epoch in training DGI')
    parser.add_argument( '--hidden', type=int, default=256, help='hidden channels in DGI') 
    parser.add_argument( '--PCA', type=int, default=1, help='run PCA or not')   
    parser.add_argument( '--cluster', type=int, default=1, help='run cluster or not')
    parser.add_argument( '--n_clusters', type=int, default=7, help='number of clusters in Kmeans, when ground truth label is not avalible.') #5 on MERFISH, 20 on Breast
    parser.add_argument( '--draw_map', type=int, default=1, help='run drawing map')
    parser.add_argument( '--diff_gene', type=int, default=0, help='Run differential gene expression analysis')
    args = parser.parse_args()

    data_fold = args.data_path+args.data_name+'/'
    generated_data_fold = args.generated_data_path + args.data_name+'/'
    if not os.path.exists(generated_data_fold):
        os.makedirs(generated_data_fold)
    adata_h5 = adata
    print(adata_h5)
    
    features = adata_preprocess(adata_h5, min_cells=args.min_cells, pca_n_comps=args.Dim_PCA)
    gene_ids = adata_h5.var['gene_ids']
    coordinates = adata_h5.obsm['spatial']

    np.save(generated_data_fold + 'features.npy', features)
    np.save(generated_data_fold + 'coordinates.npy', np.array(coordinates))


    get_adj(generated_data_fold)


    ##exocrine GCNG with normalized graph matrix 

    args.embedding_data_path = args.embedding_data_path +'/'+ args.data_name +'/'
    args.model_path = args.model_path +'/'+ args.data_name +'/'
    args.result_path = args.result_path +'/'+ args.data_name +'/'
    if not os.path.exists(args.embedding_data_path):
        os.makedirs(args.embedding_data_path) 
    if not os.path.exists(args.model_path):
        os.makedirs(args.model_path) 
    args.result_path = args.result_path+'lambdaI'+str(args.lambda_I) +'/'
    if not os.path.exists(args.result_path):
        os.makedirs(args.result_path) 
    print ('------------------------Model and Training Details--------------------------')
    print(args) 

    lambda_I = args.lambda_I
    # Parameters
    batch_size = 1  # Batch size

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    adj_0, adj, X_data = get_data(args)


    num_cell = X_data.shape[0]
    num_feature = X_data.shape[1]
    print('Adj:', adj.shape, 'Edges:', len(adj.data))
    print('X:', X_data.shape)
    
    n_clusters=7

    if args.DGI and (lambda_I>=0):
        print("-----------Deep Graph Infomax-------------")
        data_list = get_graph(adj, X_data)
        data_loader = DataLoader(data_list, batch_size=batch_size)
        DGI_model = train_DGI(args, data_loader=data_loader, in_channels=num_feature)

        for data in data_loader:
            data.to(device)
            X_embedding, _, _ = DGI_model(data)
            X_embedding = X_embedding.cpu().detach().numpy()
            X_embedding_filename =  args.embedding_data_path+'lambdaI' + str(lambda_I) + '_epoch' + str(args.num_epoch) + '_Embed_X.npy'
            np.save(X_embedding_filename, X_embedding)

    cluster_type = 'leiden' 

    print("-----------Clustering-------------")

    X_embedding_filename =  args.embedding_data_path+'lambdaI' + str(lambda_I) + '_epoch' + str(args.num_epoch) + '_Embed_X.npy'
    X_embedding = np.load(X_embedding_filename)

     
    adata = ad.AnnData(X_embedding)
    adata_spatial = sc.read_h5ad(f'{datadir}/151673.h5ad')
    adata.uns['spatial']=adata_spatial.uns['spatial']
    adata.obsm['spatial']=adata_spatial.obsm['spatial']

    gt=adata_spatial.obs['Region'].tolist()
    adata.obs['Region'] = gt
    adata.obs['Region']=adata.obs['Region'].astype('category')

    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50) # 20
    eval_resolution = res_search_fixed_clus(cluster_type, adata, n_clusters)
    
    sc.tl.leiden(adata, key_added="CCST_leiden", resolution=eval_resolution)
        
    
    obs_df = adata.obs.dropna()
    ari = adjusted_rand_score(obs_df['CCST_leiden'], obs_df['Region'])
    print('ari:{}'.format(ari))
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

    domain=adata.obs['CCST_leiden'].tolist()
    new_data.obs['CCST_{}'.format(i+1)]=domain
    new_data.obs['CCST_{}'.format(i+1)]=new_data.obs['CCST_{}'.format(i+1)].astype('category')

result_path='./'
new_data.uns['time']=during_time
new_data.uns['memory']=memory
new_data.uns['ari']=ari_list
new_data.write(f'{result_path}/CCST_151673.h5ad')


