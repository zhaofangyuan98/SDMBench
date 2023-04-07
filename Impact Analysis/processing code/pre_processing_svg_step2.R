devtools::install_github('xzhoulab/SPARK')
library('SPARK')

counts <- read.table("./spark_count_matrix.csv", check.names = F,header = TRUE, row.names = 1,sep = ",")

spatial<-read.table("./spark_spatial.csv", check.names = F,header = TRUE,sep = ",")


## extract the coordinates from the rawdata
info <- cbind.data.frame(x=spatial['x'],
                         y=spatial['y'],
                         total_counts=apply(counts,2,sum))
rownames(info) <- colnames(counts)

## filter genes and cells/spots and 
spark <- CreateSPARKObject(counts=counts, 
                           location=info[,1:2],
                           percentage = 0.1, 
                           min_total_counts = 10)

## total counts for each cell/spot
spark@lib_size <- apply(spark@counts, 2, sum)

## Take the first ten genes as an example
#spark@counts   <- spark@counts[1:10,]

## Estimating Parameter Under Null
spark <- spark.vc(spark, 
                  covariates = NULL, 
                  lib_size = spark@lib_size, 
                  num_core = 5,
                  verbose = F)

## Calculating pval
spark <- spark.test(spark, 
                    check_positive = T, 
                    verbose = F)


df_gene<-spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
result<-df_gene[order(df_gene$combined_pvalue,decreasing = TRUE),]
gene_spark<-rownames(result)[0:3000]
write.table(gene_spark,'./gene_spark.csv',row.names=FALSE)

