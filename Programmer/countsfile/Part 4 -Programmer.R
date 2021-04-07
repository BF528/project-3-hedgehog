##Installing DE based on binomial distribution  via https://bioconductor.org/packages/release/bioc/html/DESeq2.html  ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

###Reading in control csv files /project/bf528/project_3/samples/control_counts.csv   ### 
ctrl_counts<-read.csv("/projectnb/bf528/users/hedgehog/project_3/project-3-hedgehog/Programmer/countsfile/control_counts.csv")
cc_all[is.na(cc_all)] = 0
##Joining control file with treated file###

m_counts<-merge(cc_all, ctrl_counts, by= "Geneid")
countData = as.matrix(cc_all[ , -1])

coldata <- data.frame(row.names=colnames(cc_all))
condition = m_counts
dds <- DESeqDataSetFromMatrix(countData= countData, colData=coldata, )

