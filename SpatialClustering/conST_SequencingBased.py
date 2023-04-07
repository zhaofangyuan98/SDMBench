import scanpy as sc
from sklearn import metrics

import torch
import argparse
import random
import numpy as np
import pandas as pd
from src.graph_func import graph_construction
from src.utils_func import mk_dir, adata_preprocess, load_ST_file, res_search_fixed_clus, plot_clustering
from src.training import conST_training

import anndata
from sklearn import metrics
import matplotlib.pyplot as plt
import scanpy as sc
import os
import warnings
warnings.filterwarnings('ignore')

datadir = '../Data'#dataset path
memory=[1,2,3,4,5,6,7,8,9,10]
during_time=[1,2,3,4,5,6,7,8,9,10]
ari_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/151673.h5ad')


# set seed before every run
def seed_torch(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

#Run 10 repeated experiments
for j in range(10):

    import tracemalloc
    import time
 
    tracemalloc.start() 
    start_time=time.time()

    
    adata_h5 = sc.read_h5ad(f'{datadir}/151673.h5ad')
    

    parser = argparse.ArgumentParser()
    parser.add_argument('--k', type=int, default=20, help='parameter k in spatial graph')
    parser.add_argument('--knn_distanceType', type=str, default='euclidean',
                        help='graph distance type: euclidean/cosine/correlation')
    parser.add_argument('--epochs', type=int, default=200, help='Number of epochs to train.')
    parser.add_argument('--cell_feat_dim', type=int, default=32, help='Dim of PCA')
    parser.add_argument('--feat_hidden1', type=int, default=100, help='Dim of DNN hidden 1-layer.')
    parser.add_argument('--feat_hidden2', type=int, default=20, help='Dim of DNN hidden 2-layer.')
    parser.add_argument('--gcn_hidden1', type=int, default=32, help='Dim of GCN hidden 1-layer.')
    parser.add_argument('--gcn_hidden2', type=int, default=8, help='Dim of GCN hidden 2-layer.')
    parser.add_argument('--p_drop', type=float, default=0.2, help='Dropout rate.')
    parser.add_argument('--use_img', type=bool, default=False, help='Use histology images.')
    parser.add_argument('--img_w', type=float, default=0.1, help='Weight of image features.')
    parser.add_argument('--use_pretrained', type=bool, default=False, help='Use pretrained weights.')
    parser.add_argument('--using_mask', type=bool, default=False, help='Using mask for multi-dataset.')
    parser.add_argument('--feat_w', type=float, default=10, help='Weight of DNN loss.')
    parser.add_argument('--gcn_w', type=float, default=0.1, help='Weight of GCN loss.')
    parser.add_argument('--dec_kl_w', type=float, default=10, help='Weight of DEC loss.')
    parser.add_argument('--gcn_lr', type=float, default=0.01, help='Initial GNN learning rate.')
    parser.add_argument('--gcn_decay', type=float, default=0.01, help='Initial decay rate.')
    parser.add_argument('--dec_cluster_n', type=int, default=10, help='DEC cluster number.')
    parser.add_argument('--dec_interval', type=int, default=20, help='DEC interval nnumber.')
    parser.add_argument('--dec_tol', type=float, default=0.00, help='DEC tol.')

    parser.add_argument('--seed', type=int, default=0, help='random seed')
    parser.add_argument('--beta', type=float, default=100, help='beta value for l2c')
    parser.add_argument('--cont_l2l', type=float, default=0.3, help='Weight of local contrastive learning loss.')
    parser.add_argument('--cont_l2c', type=float, default= 0.1, help='Weight of context contrastive learning loss.')
    parser.add_argument('--cont_l2g', type=float, default= 0.1, help='Weight of global contrastive learning loss.')

    parser.add_argument('--edge_drop_p1', type=float, default=0.1, help='drop rate of adjacent matrix of the first view')
    parser.add_argument('--edge_drop_p2', type=float, default=0.1, help='drop rate of adjacent matrix of the second view')
    parser.add_argument('--node_drop_p1', type=float, default=0.2, help='drop rate of node features of the first view')
    parser.add_argument('--node_drop_p2', type=float, default=0.3, help='drop rate of node features of the second view')

    # ______________ Eval clustering Setting ______________
    parser.add_argument('--eval_resolution', type=int, default=1, help='Eval cluster number.')
    parser.add_argument('--eval_graph_n', type=int, default=20, help='Eval graph kN tol.') 

    params =  parser.parse_args()

    np.random.seed(params.seed)
    torch.manual_seed(params.seed)
    torch.cuda.manual_seed(params.seed)
    device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
    print('Using device: ' + device)
    params.device = device

   
    adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
    graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)
    np.save(f'./input/dlpfc/adatax.npy', adata_X)
    np.save(f'./input/dlpfc/graphdict.npy', graph_dict, allow_pickle = True)

    seed_torch(params.seed)

    data_name = 'dlpfc'
    save_root = './output/spatialLIBD/'
    data_root = '../spatialLIBD'

    params.save_path = mk_dir(f'./conST')




    adata_X = np.load(f'./input/dlpfc/adatax.npy')
    graph_dict = np.load(f'./input/dlpfc/graphdict.npy',  allow_pickle = True).item()
    params.cell_num = adata_h5.shape[0]

    n_clusters = 7
    if params.use_img:
        img_transformed = np.load('./MAE-pytorch/extracted_feature.npy')
        img_transformed = (img_transformed - img_transformed.mean()) / img_transformed.std() * adata_X.std() + adata_X.mean()
        conST_net = conST_training(adata_X, graph_dict, params, n_clusters, img_transformed)
    else:
        conST_net = conST_training(adata_X, graph_dict, params, n_clusters)
    if params.use_pretrained:
        conST_net.load_model('conST_151673.pth')
    else:
        conST_net.pretraining()
        conST_net.major_training()

    conST_embedding = conST_net.get_embedding()

    np.save(f'{params.save_path}/conST_result.npy', conST_embedding)
    # clustering
    adata_conST = anndata.AnnData(conST_embedding)
    # adata_conST.uns['spatial'] = adata_h5.uns['spatial']
    adata_conST.obsm['spatial'] = adata_h5.obsm['spatial']

    sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)

    eval_resolution = res_search_fixed_clus(adata_conST, n_clusters)
    print(eval_resolution)
    cluster_key = "conST_leiden"
    sc.tl.leiden(adata_conST, key_added=cluster_key, resolution=eval_resolution)

    domain=adata_conST.obs[cluster_key].tolist()
    adata_h5.obs['pred_{}'.format(j+1)]=domain
    adata_h5.obs['pred_{}'.format(j+1)]=adata_h5.obs['pred_{}'.format(j+1)].astype('category')
   

    obs_df = adata_h5.obs.dropna()
    ari = metrics.adjusted_rand_score(obs_df['pred_{}'.format(j+1)], obs_df['Region'])


    ari_list[j]=ari


    end_time=time.time()
    during=end_time-start_time

    size, peak = tracemalloc.get_traced_memory()

    tracemalloc.stop()
    
    memory[j]=peak /1024/1024
    during_time[j]=during
    
   
    print('memory blocks peak:{:>10.4f} MB'.format(memory[j]))
    print('time: {:.4f} s'.format(during_time[j]))
    new_data.obs['pred_{}'.format(j+1)]=adata_h5.obs['pred_{}'.format(j+1)]


result_path='./'
new_data.uns['time']=during_time
new_data.uns['memory']=memory
new_data.uns['ari']=ari_list
new_data.write(f'{result_path}/conST_151673.h5ad')


