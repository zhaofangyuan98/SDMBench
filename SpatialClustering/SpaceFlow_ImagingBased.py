import scanpy as sc
from sklearn import metrics
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import squidpy as sq
import SpaceFlow
import scanpy as sc
import pandas as pd

from SpaceFlow import SpaceFlow 

datadir = '../Data'#dataset path
memory=[1,2,3,4,5,6,7,8,9,10]
during_time=[1,2,3,4,5,6,7,8,9,10]
ari_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/osmfish.h5ad')


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

    
    adata = sc.read_h5ad(f'{datadir}/osmfish.h5ad')
    print(adata)

    adata.var_names_make_unique()
    sf = SpaceFlow.SpaceFlow(adata=adata)

    #preprocess
    sf.preprocessing_data(n_top_genes=3000)

    sf.train(spatial_regularization_strength=0.1, 
         z_dim=50, 
         lr=1e-3, 
         epochs=1000, 
         max_patience=50, 
         min_stop=100, 
         random_seed=42, 
         gpu=1, 
         regularization_acceleration=True, 
         edge_subset_sz=1000000)
    
    n_clusters=11
    sc.pp.neighbors(adata, n_neighbors=50)
    eval_resolution = res_search_fixed_clus(adata, n_clusters)

    sf.segmentation(domain_label_save_filepath="./osmfish_domains_{}.csv".format(i+1), 
                n_neighbors=50, 
                resolution=eval_resolution)

    sf.plot_segmentation(segmentation_figure_save_filepath="./osmfish_domain_segmentation_{}.pdf".format(i+1), 
                     colormap="tab20", 
                     scatter_sz=1., 
                     rsz=4., 
                     csz=4., 
                     wspace=.4, 
                     hspace=.5, 
                     left=0.125, 
                     right=0.9, 
                     bottom=0.1, 
                     top=0.9)
    
    pred=pd.read_csv('./osmfish_domains_{}.csv'.format(i+1),header=None)
    pred_list=pred.iloc[:,0].to_list()
    
    new_data.obs['pred_{}'.format(i+1)] = np.array(pred_list)

    obs_df = new_data.obs.dropna()
    ari = metrics.adjusted_rand_score(obs_df['Region'], obs_df['pred_{}'.format(i+1)])
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


result_path='./'
new_data.uns['time']=during_time
new_data.uns['memory']=memory
new_data.uns['ari']=ari_list
new_data.write(f'{result_path}/SpaceFlow_osmfish.h5ad')


