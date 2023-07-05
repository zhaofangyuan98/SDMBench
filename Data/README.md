# Datasets

There are some spatial transcriptome datasets used in our benchmark:

***osmfish.h5ad*** is a spatial transcriptomics dataset of mouse somatosensory cortex. obs['Region'] is the cell type annotation of each spot, and we use this info as the ground truth. We have removed the "Excluded" region of this dataset.

***bass_count_matrix_osmfish.csv*** is the gene expression matrix file used in the osmfish dataset for the BASS method.

***bass_spatial_osmfish.csv*** is the spatial location file used in the osmfish dataset for the BASS method.

Our reprocessed versions of all datasets are publicly available as h5ad format on Figshare (<https://figshare.com/projects/SDMBench/163942>) or website (<http://sdmbench.drai.cn/>).
