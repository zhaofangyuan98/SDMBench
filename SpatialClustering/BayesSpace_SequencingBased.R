rm(list = ls())

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))
BiocManager::install("BayesSpace")

devtools::install_github("cellgeni/sceasy")

install.packages('reticulate')
install.packages('Seurat')
library(reticulate)
library(sceasy)
library(Seurat)

#replace your python path
Sys.setenv(RETICULATE_PYTHON="D:\\anaconda3\\python.exe")
use_python("D:\\anaconda3\\python.exe")
py_config()

sceasy::convertFormat("../Data/151673.h5ad", from="anndata", to="seurat",
                      outFile='../Data/151673.rds')

library(scater)
library(SeuratData)
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
library(patchwork)
library(BayesSpace)
library(ggplot2)

runningtime<-c(1,2,3,4,5,6,7,8,9,10)

#Run 10 repeated experiments
for (i in 1:10) {
  
  gc(reset = T)
  
  timestart<-Sys.time()
  
  mbsp_temp<-readRDS(file = '../Data/151673.rds')
  mbsp <-as.SingleCellExperiment(mbsp_temp)
 
  set.seed(101)
  dec <- scran::modelGeneVar(mbsp)
  top <- scran::getTopHVGs(dec, n = 2000)
  
  set.seed(102)
  mbsp <- scater::runPCA(mbsp, subset_row=top)
  
  ## Add BayesSpace metadata
  mbsp <- spatialPreprocess(mbsp, platform="Visium", skip.PCA=TRUE)
  
  
  #Clustering with BayesSpace
  q <- 7  # Number of clusters[5,6,7,8,9,10,11,12,13,14]
  d <- 15  # Number of PCs
  
  ## Run BayesSpace clustering
  set.seed(104)
  mbsp <- spatialCluster(mbsp, q=q, d=d, platform='Visium',
                          nrep=50000, gamma=3, save.chain=TRUE)
  
  labels <- mbsp$spatial.cluster
  label_table <-data.frame(labels)
  write.table(label_table,paste('./BayesSpace_dataset_',i,'.csv'),row.names=FALSE)

  timeend<-Sys.time()
  runningtime[i] <-timeend-timestart
  
  print('runningtime:')
  print(runningtime[i])
  
  #内存
  gc()
  memInfo1 <- gc()
  memInfo1[11] 
  memInfo1[12] 
  
  gc(reset=TRUE)
  memInfo2 <- gc()
  memInfo2[11] 
  memInfo2[12] 
  peak<-memInfo1[12]-memInfo2[12]
  print(peak)
  
}








