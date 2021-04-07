#Mackenzie Knox
#BF528 Project 3
#Biologist

#Imports
library(readr)
library(dplyr)

#set working directory for reading csv
setwd("/projectnb/bf528/users/hedgehog/project_3/programmer/csv_files")

#read for filtering and heatmapping
ahr <- read.csv("CMC_AHR_deseq_norm_counts.csv")
cxr<- read.csv("CORN__CARPXR_deseq_norm_counts.csv")
ctc <- read.csv("SALINE_Cytotoxic_deseq_norm_counts.csv")

setwd("/projectnb/bf528/users/hedgehog/project_3/biologist") #just in case

colnames(ahr)[2:7] <- paste(colnames(ahr)[2:7], "AhR")
colnames(ctc)[2:7] <- paste(colnames(ctc)[2:7], "Cytotoxic")
colnames(cxr)[2:7] <- paste(colnames(cxr)[2:7], "CARPXR")


#Heatmapping

#funciton for cov filter
cv <- function(row){sd(row)/mean(row)}
cutoff <- 0.182

#Merge dataframes for heatmap

df <- merge(ahr, ctc, by="X", all = TRUE)
df <- merge(df, cxr, by="X", all=TRUE)
rownames(df) <- df$X
df <- df[,-1]


#filtering
dfm <- median(apply(df, 1,FUN=mean, na.rm=T))
dfmean <- df[apply(df, 1,FUN=mean, na.rm=T)<dfm,]


dfcov <- df[apply(df, 1, FUN=cv)<cutoff,]
dfmcov <- (apply(dfcov, 1,FUN=mean, na.rm=T))
dfmeancov <- dfcov[apply(dfcov, 1,FUN=mean, na.rm=T)<dfmcov[2],]

dfmcov <- apply(dfcov, 1,FUN=mean, na.rm=T)
#Median is third
dfmeanm <- dfcov[apply(dfcov, 1,FUN=mean, na.rm=T)<dfmcov[3],]

#Create Heatmap
heatmap(data.matrix(dfmeanm), margins=c(10,0))

library(gplots)
heatmap.2(data.matrix(dfmeanm), 
          main='Gene Expression',
          trace='none', density.info = 'none',
          key.xlab='Gene Expression Level', scale='row', margins=c(10,3))


setwd('/projectnb/bf528/users/hedgehog/project_3/programmer/csv_files')

# subset the data to get the significant DE genes (adjusted p-value <0.05) 
AHR<- read.csv("resCMC_AHR_deseq_results.csv", header=TRUE)
CARPXR<-read.csv("resCORN_CARPXR_deseq_results.csv", header=TRUE)
CYTO<-read.csv("resSALINE_Cytotoxic_deseq_results.csv", header=TRUE)
sig_AHR <- subset(AHR, pvalue < 0.05)
sig_CARPXR <- subset(CARPXR, pvalue < 0.05)
sig_CYTO <- subset(CYTO, pvalue < 0.05)
head(sig_AHR)

#subset log2foldchange > 1.5
sig_AHR1 <- subset(sig_AHR, abs(log2FoldChange) > 1.5)
sig_CARPXR1 <- subset(sig_CARPXR, abs(log2FoldChange) > 1.5)
sig_CYTO1 <- subset(sig_CYTO, abs(log2FoldChange) > 1.5)
#dim(sig_AHR1)
#dim(sig_CARPXR1)
#dim(sig_CYTO1)

#write to csv files
sig_AHR1<- cbind(sig_AHR1)
sig_CARPXR1 <- cbind(sig_CARPXR1)
sig_CYTO1<- cbind(sig_CYTO1)

setwd('/projectnb/bf528/users/hedgehog/project_3/biologist')
write.csv(sig_AHR1, "sig_AHR.csv")
write.csv(sig_CARPXR1, "sig_CARPXR.csv")
write.csv(sig_CYTO1, "sig_Cytotoxic.csv")
