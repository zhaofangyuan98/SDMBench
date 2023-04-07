rm(list = ls())

if(!require(devtools))
  install.packages(devtools)
devtools::install_github("zhengli09/BASS")

library(BASS)
library(Matrix)

counts <- read.table("../Data/bass_count_matrix_osmfish.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts <- as(as.matrix(counts), "dgCMatrix")

spatial<-read.csv("../Data/bass_spatial_osmfish.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")

runningtime<-c(1,2,3,4,5,6,7,8,9,10)

#Run 10 repeated experiments
for (i in 1:10) {
  
  gc(reset = T)
  
  timestart<-Sys.time()
  
  
  smp <- "osmfish"
  
  # hyper-parameters
  # We set the number of cell types to a relatively large
  # number (20) to capture the expression heterogeneity.
  C <- 20
  # number of spatial domains
  R <- 11
  
  
  set.seed(0)
  # Set up BASS object
  BASS <- createBASSObject(counts, spatial, C = C, R = R,
                           beta_method = "SW")
  
  # Data pre-processing:
  # 1.Library size normalization followed with a log2 transformation
  # 2.Select top 3000 spatially expressed genes with SPARK-X
  # 3.Dimension reduction with PCA
  BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE,
                          geneSelect = "sparkx", nSE = 3000, doPCA = TRUE, 
                          scaleFeature = FALSE, nPC = 20)
  
  # Run BASS algorithm
  BASS <- BASS.run(BASS)
  
  # post-process posterior samples:
  # 1.Adjust for label switching with the ECR-1 algorithm
  # 2.Summarize the posterior samples to obtain the spatial domain labels
  BASS <- BASS.postprocess(BASS)
  
  zlabels <- BASS@results$z # spatial domain labels
  
  label_table <-data.frame(zlabels)
  write.table(label_table,paste('./domain_',i,'.csv'),row.names=FALSE)
  
  timeend<-Sys.time()
  runningtime[i] <-timeend-timestart
  
  print('runningtime:')
  print(runningtime[i])
  
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








