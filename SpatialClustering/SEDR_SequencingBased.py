import scanpy as sc
from sklearn import metrics

import torch
import argparse
import warnings
import numpy as np
import pandas as pd
from src.graph_func import graph_construction
from src.utils_func import mk_dir, adata_preprocess, load_ST_file
import anndata
from src.SEDR_train import SEDR_Train
from sklearn import metrics
import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings('ignore')

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

    adata = sc.read_h5ad(f'{datadir}/151673.h5ad')
    
    torch.cuda.cudnn_enabled = False
    np.random.seed(0)
    torch.manual_seed(0)
    torch.cuda.manual_seed(0)
    device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
    print('===== Using device: ' + device)

    # ################ Parameter setting
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', type=int, default=10, help='parameter k in spatial graph')
    parser.add_argument('--knn_distanceType', type=str, default='euclidean',
                        help='graph distance type: euclidean/cosine/correlation')
    parser.add_argument('--epochs', type=int, default=200, help='Number of epochs to train.')
    parser.add_argument('--cell_feat_dim', type=int, default=32, help='Dim of PCA')#151673的200变为32
    parser.add_argument('--feat_hidden1', type=int, default=100, help='Dim of DNN hidden 1-layer.')
    parser.add_argument('--feat_hidden2', type=int, default=20, help='Dim of DNN hidden 2-layer.')
    parser.add_argument('--gcn_hidden1', type=int, default=32, help='Dim of GCN hidden 1-layer.')
    parser.add_argument('--gcn_hidden2', type=int, default=8, help='Dim of GCN hidden 2-layer.')
    parser.add_argument('--p_drop', type=float, default=0.2, help='Dropout rate.')
    parser.add_argument('--using_dec', type=bool, default=True, help='Using DEC loss.')
    parser.add_argument('--using_mask', type=bool, default=False, help='Using mask for multi-dataset.')
    parser.add_argument('--feat_w', type=float, default=10, help='Weight of DNN loss.')
    parser.add_argument('--gcn_w', type=float, default=0.1, help='Weight of GCN loss.')
    parser.add_argument('--dec_kl_w', type=float, default=10, help='Weight of DEC loss.')
    parser.add_argument('--gcn_lr', type=float, default=0.01, help='Initial GNN learning rate.')
    parser.add_argument('--gcn_decay', type=float, default=0.01, help='Initial decay rate.')
    parser.add_argument('--dec_cluster_n', type=int, default=10, help='DEC cluster number.')
    parser.add_argument('--dec_interval', type=int, default=20, help='DEC interval nnumber.')
    parser.add_argument('--dec_tol', type=float, default=0.00, help='DEC tol.')
    # ______________ Eval clustering Setting ______________
    parser.add_argument('--eval_resolution', type=int, default=1, help='Eval cluster number.')
    parser.add_argument('--eval_graph_n', type=int, default=20, help='Eval graph kN tol.') 
    parser.add_argument( '--name', type=str)
    parser.add_argument( '--n_cluster', type=int)
    parser.add_argument( '--input_dir', type=str)
    parser.add_argument( '--output_dir', type=str)

    params = parser.parse_args()
    params.device = device

    
    n_clusters = 7


    # ################## Load data
    # adata_h5 = load_ST_file(file_fold=file_fold)
    adata_h5=adata
    adata_h5.var_names_make_unique()

    adata_X = adata_preprocess(adata_h5, min_cells=5, pca_n_comps=params.cell_feat_dim)
    graph_dict = graph_construction(adata_h5.obsm['spatial'], adata_h5.shape[0], params)
    params.save_path = mk_dir(f'./')

    params.cell_num = adata_h5.shape[0]
    print('==== Graph Construction Finished')

    # ################## Model training
    sedr_net = SEDR_Train(adata_X, graph_dict, params)
    if params.using_dec:
        sedr_net.train_with_dec()
    else:
        sedr_net.train_without_dec()
    sedr_feat, _, _, _ = sedr_net.process()

    np.savez(f'./SEDR_result.npz', sedr_feat=sedr_feat, params=params)

    # ################## Result plot
    adata_sedr = anndata.AnnData(sedr_feat)
    
    gt=adata_h5.obs['Region'].tolist()
    adata_sedr.obs['Region'] = gt
    adata_sedr.obs['Region']=adata_sedr.obs['Region'].astype('category')

    adata_sedr.obsm['spatial'] = adata_h5.obsm['spatial']
    sc.pp.neighbors(adata_sedr, n_neighbors=params.eval_graph_n)
    eval_resolution = res_search_fixed_clus(adata_sedr, n_clusters)
    sc.tl.leiden(adata_sedr, key_added="SEDR_leiden", resolution=eval_resolution)
    
    obs_df = adata_sedr.obs.dropna()
    print(obs_df)
    ari = metrics.adjusted_rand_score(obs_df['Region'], obs_df['SEDR_leiden'])

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
    
    domain=adata_sedr.obs['SEDR_leiden'].tolist()
    new_data.obs['pred_{}'.format(i+1)]=domain
    new_data.obs['pred_{}'.format(i+1)]=new_data.obs['pred_{}'.format(i+1)].astype('category')
    
result_path='./'
new_data.uns['time']=during_time
new_data.uns['memory']=memory
new_data.uns['ari']=ari_list
new_data.write(f'{result_path}/SEDR_151673.h5ad')


