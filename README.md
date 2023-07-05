# SDMBench
We collected **87** real spatial transcriptome data and simulated **175** data to benchmark **13** methods for identifying spatial domains. **13** methods are included **Louvain, Leiden, SpaGCN, BayesSpace, StLearn, SEDR, CCST, SCAN-IT, STAGATE, SpaceFlow, conST, BASS, DeepST**.

## Overview
The main work is as follows.

1)**13** computational methods were benchmarked on **34** real data form different spatial transcriptome technologies (10x Visium, Stero-Seq, BaristaSeq, MERFISH, osmFISH, STARmap, STARmap*). The benchmark encompassed 10 different metrics to assess the methods' performance in terms of accuracy, spatial continuity, marker gene detection, scalability, and robustness. The metrics included NMI, HOM, COM, CHAOS, PAS, ASW, Moran's I, Geary's C, time, and memory.
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




| :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: |
| NMI | HOM | COM | CHAOS | PAS | ASW | Moran's I | Geary's C | time | memory |


2)Furthermore, with **175** simulated datasets, we examined the robustness of these 13 methods against four different factors (different gene expression matrix sparsity, spatial resolutions and gene numbers, levels of noise), and assessed the impact of pre- and post-processing steps on performance. 

3)We identified limitations within current methodologies. These limitations became apparent during the testing of various spatial clustering methods on additional **22** data containing small and non-continuous tissue domains, and during multi-slice analysis on another large-scale dataset containing **31** data slices. We proposed a "divide and conquer" strategy to make the latter challenging task effectively solvable.

4)We provided user guidance for selecting the optimal spatial clustering methods based on data characteristics(See Paper Section ”Overall performance” and Fig. 4, 5). [The tutorial](https://github.com/zhaofangyuan98/SDMBench/tree/main/Tutorial) serves as an illustrative example that demonstrates the process of benchmarking new methods against existing ones. Additionally, we have developed [a website interface](http://sdmbench.drai.cn/) that facilitates the benchmarking of new methods by comparing them with established approaches.

To make it easier for users to reproduce these computational methods, we plan to create a separate docker on [Docker](https://github.com/zhaofangyuan98/SDMBench/tree/main/Docker) for each method. We are actively updating this repository! 
