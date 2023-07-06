# SDMBench
We collected **87** real spatial transcriptome data and simulated **175** data to benchmark **13** methods for identifying spatial domains. 
 
**13** methods are included:
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
- **DeepST**:《DeepST: A versatile graph contrastive learning framework for spatially informed clustering, integration, and deconvolution of spatial transcriptomics》

## Overview
The main work is as follows.

### 1. Benchmark on different spatial transcriptome data
**13** computational methods were benchmarked on **34** real data form different spatial transcriptome technologies (10x Visium, Stero-Seq, BaristaSeq, MERFISH, osmFISH, STARmap, STARmap*). The benchmark encompassed 10 different metrics to assess the methods' performance in terms of accuracy, spatial continuity, marker gene detection, scalability, and robustness. The metrics included NMI, HOM, COM, CHAOS, PAS, ASW, Moran's I, Geary's C, time, and memory.

<table>
    <tr>
      <td colspan="3">Accuracy</td>   
      <td colspan="3">Continuity</td>
      <td colspan="2">Maker score</td>
      <td colspan="2">Scalability</td>
    </tr>
    <tr>
      <td>NMI</td>
      <td>HOM</td>
      <td>COM</td>
      <td>CHAOS</td>
      <td>PAS</td>
      <td>ASW</td>
      <td>Moran's I</td>
      <td>Geary's C</td>
      <td>time</td>
      <td>memory</td>
    </tr>
</table>

### 2. The robustness against four different factors
Furthermore, with **175** simulated datasets, we examined the robustness of these 13 methods against four different factors ***(different gene expression matrix sparsity, spatial resolutions and gene numbers, levels of noise)***, and assessed the impact of ***pre- and post-processing steps*** on performance. 

### 3. Limitations within current methods
We identified limitations within current methodologies. These limitations became apparent during the testing of various spatial clustering methods on additional **22** data containing small and non-continuous tissue domains, and during multi-slice analysis on another large-scale dataset containing **31** data slices. We proposed ***a "divide and conquer" strategy*** to make the latter challenging task effectively solvable.

### 4. Guidance and tutorial
We provided ***user guidance*** for selecting the optimal spatial clustering methods (See Paper Section ”Overall performance” and Fig. 4, 5). [***The tutorial***](https://github.com/zhaofangyuan98/SDMBench/tree/main/Tutorial) serves as an illustrative example that demonstrates the process of benchmarking new methods against existing ones. Additionally, we have developed [***a website interface***](http://sdmbench.drai.cn/) that facilitates the benchmarking of new methods by comparing them with established approaches.

To make it easier for users to reproduce these computational methods, we plan to create a separate docker on [Docker](https://github.com/zhaofangyuan98/SDMBench/tree/main/Docker) for each method. We are actively updating this repository! 

## Tutorial

If you want to analysis your own data, [***The tutorial***](https://github.com/zhaofangyuan98/SDMBench/tree/main/Tutorial) is an example demonstrating how to compare performance of new methods against existing ones. You can run the jupyter notebook of [tutorial.ipynb](https://github.com/zhaofangyuan98/SDMBench/blob/main/Tutorial/tutorial.ipynb) to reproduce it, and [new_method.txt](https://github.com/zhaofangyuan98/SDMBench/blob/main/Tutorial/new_method.txt) provides a format about your results file. Before you can run tutorial.ipynb, you need to install the SDMBench package (see the directory [***SDMBench***](https://github.com/zhaofangyuan98/SDMBench/tree/main/SDMBench)).

We also created a Spatial DoMain Benchmark (SDMBench) website, accessible at [http://sdmbench.drai.cn/], where you can easily benchmark your new methods against existing ones. We provided a tutorial and video (see the [“tutorial”](http://sdmbench.drai.cn/Tutorial/) in the website) for how to use the SDMBench website with examples.

## Datasets

Our reprocessed versions of all datasets are publicly available as h5ad format, for convenience datasets can be downloaded from <http://sdmbench.drai.cn/> and on Figshare (<https://figshare.com/projects/SDMBench/163942>)

For more details, please see our paper, thank you very much.
