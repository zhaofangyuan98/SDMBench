rm(list = ls())

library(BASS)
library(Seurat)
library(tidyverse)


library(Matrix)


counts_0 <- read.table("F:/BASS/aging/bass_count_matrix_MsBrainAgingSpatialDonor_12_1.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts_0 <- as(as.matrix(counts_0), "dgCMatrix")

counts_1 <- read.table("F:/BASS/aging/bass_count_matrix_MsBrainAgingSpatialDonor_12_0.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts_1 <- as(as.matrix(counts_1), "dgCMatrix")

counts_2 <- read.table("F:/BASS/aging/bass_count_matrix_MsBrainAgingSpatialDonor_11_2.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts_2 <- as(as.matrix(counts_2), "dgCMatrix")

counts_3 <- read.table("F:/BASS/aging/bass_count_matrix_MsBrainAgingSpatialDonor_11_1.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts_3 <- as(as.matrix(counts_3), "dgCMatrix")

counts_4 <- read.table("F:/BASS/aging/bass_count_matrix_MsBrainAgingSpatialDonor_11_0.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts_4 <- as(as.matrix(counts_4), "dgCMatrix")

counts_5 <- read.table("F:/BASS/aging/bass_count_matrix_MsBrainAgingSpatialDonor_10_2.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts_5 <- as(as.matrix(counts_5), "dgCMatrix")

counts_6 <- read.table("F:/BASS/aging/bass_count_matrix_MsBrainAgingSpatialDonor_10_1.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts_6 <- as(as.matrix(counts_6), "dgCMatrix")

counts_7 <- read.table("F:/BASS/aging/bass_count_matrix_MsBrainAgingSpatialDonor_10_0.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts_7 <- as(as.matrix(counts_7), "dgCMatrix")

counts_8 <- read.table("F:/BASS/aging/bass_count_matrix_MsBrainAgingSpatialDonor_9_2.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts_8 <- as(as.matrix(counts_8), "dgCMatrix")

counts_9 <- read.table("F:/BASS/aging/bass_count_matrix_MsBrainAgingSpatialDonor_9_1.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")
counts_9 <- as(as.matrix(counts_9), "dgCMatrix")

counts <- list()
counts[['000']] <- counts_0
counts[['111']] <- counts_1
counts[['222']] <- counts_2
counts[['333']] <- counts_3
counts[['444']] <- counts_4
counts[['555']] <- counts_5
counts[['666']] <- counts_6
counts[['777']] <- counts_7
counts[['888']] <- counts_8
counts[['999']] <- counts_9

spatial_0<-read.csv("F:/BASS/aging/bass_spatial_MsBrainAgingSpatialDonor_12_1.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")
spatial_1<-read.csv("F:/BASS/aging/bass_spatial_MsBrainAgingSpatialDonor_12_0.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")
spatial_2<-read.table("F:/BASS/aging/bass_spatial_MsBrainAgingSpatialDonor_11_2.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")
spatial_3<-read.table("F:/BASS/aging/bass_spatial_MsBrainAgingSpatialDonor_11_1.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")
spatial_4<-read.csv("F:/BASS/aging/bass_spatial_MsBrainAgingSpatialDonor_11_0.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")
spatial_5<-read.csv("F:/BASS/aging/bass_spatial_MsBrainAgingSpatialDonor_10_2.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")
spatial_6<-read.csv("F:/BASS/aging/bass_spatial_MsBrainAgingSpatialDonor_10_1.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")
spatial_7<-read.csv("F:/BASS/aging/bass_spatial_MsBrainAgingSpatialDonor_10_0.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")
spatial_8<-read.csv("F:/BASS/aging/bass_spatial_MsBrainAgingSpatialDonor_9_2.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")
spatial_9<-read.csv("F:/BASS/aging/bass_spatial_MsBrainAgingSpatialDonor_9_1.csv", check.names = F,header = TRUE,row.names = 1,sep = ",")


spatial<- list()
spatial[['000']] <- spatial_0
spatial[['111']] <- spatial_1
spatial[['222']] <- spatial_2
spatial[['333']] <- spatial_3
spatial[['444']] <- spatial_4
spatial[['555']] <- spatial_5
spatial[['666']] <- spatial_6
spatial[['777']] <- spatial_7
spatial[['888']] <- spatial_8
spatial[['999']] <- spatial_9




timestart<-Sys.time()

# starmap_cnts, starmap_info
#load("F:/BASS/MERFISH_Animal1.RData")

#Load data and set hyper-parameters
smps <- c("000","111","222","333","444","555","666","777","888","999")
cnts <- counts[smps] # a list of gene expression count matrices
xys <- lapply(spatial[smps], function(info.i){
  info.i$x <- info.i$x - min(info.i$x)
  info.i$y <- info.i$y - min(info.i$y)
  as.matrix(info.i[, c("x", "y")])
}) # a list of spatial coordinates matrices
# hyper-parameters
C <- 20 # number of cell types
R <- 8 # number of spatial domains


#Run BASS
set.seed(0)
# Set up BASS object
BASS <- createBASSObject(cnts, xys, C = C, R = R, beta_method = "SW")

# Data pre-processing:
# 1.Library size normalization followed with a log2 transformation
# 2.Dimension reduction with PCA after standardizing all the genes
# 3.Batch effect adjustment using the Harmony package
BASS <- BASS.preprocess(BASS, doLogNormalize = FALSE,doBatchCorrect = FALSE,
                        doPCA = TRUE, scaleFeature = FALSE, nPC = 20)

# Run BASS algorithm
BASS <- BASS.run(BASS)

# The spatial parameter beta has converged
# after checking the trace plot of beta
plot(1:BASS@burnin, BASS@samples$beta, xlab = "Iteration", 
     ylab = expression(beta), type = "l")

# post-process posterior samples:
# 1.Adjust for label switching with the ECR-1 algorithm
# 2.Summarize the posterior samples to obtain the cell type labels, spatial 
#   domain labels, and cell type proportion matrix estimate
BASS <- BASS.postprocess(BASS)

clabels <- BASS@results$c # cell type clusters
zlabels <- BASS@results$z # spatial domain labels


#I added!
write.table(zlabels[0],file='F:/BASS/domain_12_1.csv',row.names=FALSE)
write.table(zlabels[1],file='F:/BASS/domain_12_0.csv',row.names=FALSE)
write.table(zlabels[2],file='F:/BASS/domain_11_2.csv',row.names=FALSE)
write.table(zlabels[3],file='F:/BASS/domain_11_1.csv',row.names=FALSE)
write.table(zlabels[4],file='F:/BASS/domain_11_0.csv',row.names=FALSE)
write.table(zlabels[5],file='F:/BASS/domain_10_2.csv',row.names=FALSE)
write.table(zlabels[6],file='F:/BASS/domain_10_1.csv',row.names=FALSE)
write.table(zlabels[7],file='F:/BASS/domain_10_0.csv',row.names=FALSE)
write.table(zlabels[8],file='F:/BASS/domain_9_2.csv',row.names=FALSE)
write.table(zlabels[9],file='F:/BASS/domain_9_1.csv',row.names=FALSE)


timeend<-Sys.time()
runningtime <-timeend-timestart

print('runningtime:')
print(runningtime)