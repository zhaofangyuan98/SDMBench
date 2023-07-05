## SpatialClustering

There are a total of ***13*** spatial clustering methods in our study, each with its unique approach to spatial analysis. The methods included are Louvain, Leiden, SpaGCN (with and without histological image, denoted as SpaGCN(HE) and SpaGCN respectively), BayesSpace, StLearn, SEDR, CCST, SCAN-IT, STAGATE, SpaceFlow, conST, BASS, and DeepST.

- **Louvain**:《SCANPY: large-scale single-cell gene expression data analysis》
- **Leiden**:《From Louvain to Leiden: guaranteeing well-connected communities》
- **SpaGCN**: 《SpaGCN: Integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network》
- **BayesSpace**:《Spatial transcriptomics at subspot resolution with BayesSpace》
- **StLearn**: 《stLearn: integrating spatial location, tissue morphology and gene expression to find cell types, cell-cell interactions and spatial trajectories within undissociated tissues》
- **SEDR**:《Unsupervised Spatially Embedded Deep Representation of Spatial Transcriptomics》
- **CCST**:《CCST: Cell clustering for spatial transcriptomics data with graph neural network》
- **SCAN-IT**:《SCAN-IT: Domain segmentation of spatial transcriptomics images by graph neural network》
- **STAGATE**:《Deciphering spatial domains from spatially resolved transcriptomics with an adaptive graph attention auto-encoder》
- **SpaceFlow**:《Identifying multicellular spatiotemporal organization of cells with SpaceFlow》
- **conST**:《conST: an interpretable multi-modal contrastive learning framework for spatial transcriptomics》
- **BASS**: 《BASS: multi‑scale and multi‑sample analysis enables accurate cell type clustering and spatial domain detection in spatial transcriptomic studies》
- **DeepST**:《Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST》

Each evaluation in our study involved ***10*** replicated runs of each method. We performed spatial clustering analysis on two datasets as representative examples: Data 9, which is based on sequencing using the 10x Visium technology, specifically the DLPFC dataset of sample 151673, and Data 21, which is based on imaging using the osmFISH technology. 
 
The data can be download from <http://sdmbench.drai.cn/>.
