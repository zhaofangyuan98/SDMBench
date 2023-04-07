import scanpy as sc
from sklearn import metrics

datadir = '../Data'#dataset path
memory=[1,2,3,4,5,6,7,8,9,10]
during_time=[1,2,3,4,5,6,7,8,9,10]
ari_list=[1,2,3,4,5,6,7,8,9,10]

new_data = sc.read_h5ad(f'{datadir}/151673.h5ad')

for i in range(10):

    import tracemalloc
    import time
 
    tracemalloc.start()
    start_time=time.time()

    adata = sc.read_h5ad(f'{datadir}/151673.h5ad')
    
    import stlearn as st
    import scanpy as sc
    st.settings.set_figure_params(dpi=180)

    n_cluster=7

    # load data
    data = st.Read10X("./spatialLIBD/151673/")
    adata.obs['imagecol']=data.obs['imagecol']
    adata.obs['imagerow']=data.obs['imagerow']
    adata.uns['spatial']=data.uns['spatial']


    # pre-processing for gene count table
    st.pp.filter_genes(adata,min_cells=1)
    st.pp.normalize_total(adata)
    st.pp.log1p(adata)

    # run PCA for gene expression data
    st.em.run_pca(adata,n_comps=15)



    if not os.path.isdir("./151673/preprocessing/"):
        os.mkdir("./151673/preprocessing/")

    # # pre-processing for spot image
    st.pp.tiling(adata, "./151673/preprocessing/")


    # this step uses deep learning model to extract high-level features from tile images
    # may need few minutes to be completed
    st.pp.extract_feature(adata)

    # stSME
    st.spatial.SME.SME_normalize(adata, use_data="raw", weights="physical_distance")
    adata_ = adata.copy()
    adata_.X = adata_.obsm['raw_SME_normalized']

    st.pp.scale(adata_)
    st.em.run_pca(adata_,n_comps=15)

    # K-means clustering on stSME normalised PCA
    import matplotlib.pyplot as plt
    st.tl.clustering.kmeans(adata_, n_clusters=n_cluster, use_data="X_pca", key_added="X_pca_kmeans")

    gt=adata.obs['Region'].tolist()
    adata_.obs['Region'] = gt
    adata_.obs['Region']=adata_.obs['Region'].astype('category')

    obs_df = adata_.obs.dropna()
    ari = metrics.adjusted_rand_score(obs_df['Region'], obs_df['X_pca_kmeans'])
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

    new_data.obs['stLearn']=adata_.obs['X_pca_kmeans']

result_path='./'
new_data.uns['time']=during_time
new_data.uns['memory']=memory
new_data.uns['ari']=ari_list
new_data.write(f'{result_path}/stLearn_151673.h5ad')


