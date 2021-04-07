library(tidyverse)
library(dplyr)
library(knitr)
library(DescTools)
library(purrr)
library(data.table)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
library(DESeq2)
library(apeglm)
library(RColorBrewer)
library(tibble)
library(ggplot2)


###Change working directory
setwd("/projectnb2/bf528/users/hedgehog/project_3/project-3-hedgehog/Programmer/countsfile/counts_only")


#Merging all .txt files

SRR966<-read.table("SRR1177966.txt", header = TRUE)%>%select(1,7)
SRR969<-read.table("SRR1177969.txt", header = TRUE)%>%select(1,7)
SRR970<-read.table("SRR1177970.txt", header = TRUE)%>%select(1,7)
SRR993<-read.table("SRR1177993.txt", header = TRUE)%>%select(1,7)
SRR994<-read.table("SRR1177994.txt", header = TRUE)%>%select(1,7)
SRR995<-read.table("SRR1177995.txt", header = TRUE)%>%select(1,7)
SRR998<-read.table("SRR1177998.txt", header = TRUE)%>%select(1,7)
SRR8001<-read.table("SRR1178001.txt", header = TRUE)%>%select(1,7)
SRR8003<-read.table("SRR1178003.txt", header = TRUE)%>%select(1,7)

#renaming count columns with file names
names(SRR966)[2] <- "SRR1177966"
names(SRR969)[2] <- "SRR1177969"
names(SRR970)[2] <-"SRR1177970"
names(SRR993)[2] <-"SRR1177993"
names(SRR994)[2] <-"SRR1177994"
names(SRR995)[2] <-"SRR1177995"
names(SRR998)[2] <-"SRR1177998"
names(SRR8001)[2] <-"SRR1178001"
names(SRR8003)[2] <-"SRR1178003"

#use MultMerge to join greater than 2 dataframes 

merg_cnts <- MultMerge(SRR966,SRR969,SRR970,SRR993,SRR994,SRR995,SRR998,SRR8001,SRR8003, by = "Geneid")

#create csv files

write.csv(merg_cnts, "merg_cnts.csv",row.names = FALSE)


# 3.5  Creating box plots for merged samples (need to finish this one to include log scale)
## Boxplot added for 3.5
merg_cnts2 <- subset(merg_cnts, rowSums(merg_cnts== 0) == 0)
rownames(merg_cnts2)<-merg_cnts2$Geneid
merg_cnts2<-merg_cnts2[,-1]
as.matrix(merg_cnts$Geneid)

xnam = colnames(merg_cnts2)

bx_cnts<-boxplot(merg_cnts2, main = "Merged Count Data", 
                 xlab = "Runs",   lex.order = TRUE, sep = "", notch = FALSE,
                 cex.names=0.2,cex.axis=0.6,
                 ylab = "Counts", col = rainbow(9), log = "y", range = 2)
which.max(bx_cnts$out)
which.min(bx_cnts$out)
min<-bx_cnts$out[728]


#Part 4

# reading  in metadata files
meta1<-read.csv("/projectnb/bf528/users/hedgehog/project_3/project-3-hedgehog/Programmer/countsfile/toxgroup_2_rna_info.csv")
rownames(meta1) <- meta1$Run
meta1<-meta1[,-1]

#reading in control counts file
ctrl_cnts<-read.csv("/projectnb/bf528/users/hedgehog/project_3/project-3-hedgehog/Programmer/countsfile/control_counts.csv")

##View control rows with this code: colnames(ctrl_cnts)
##View meta columns with this code: rownames(meta1)

#Join counts with controls using Geneid
cc_ctrl<- full_join(merg_cnts, ctrl_cnts, by="Geneid")
cc_ctrl <- subset(cc_ctrl, rowSums(cc_ctrl== 0) == 0)
rownames(cc_ctrl)<-cc_ctrl$Geneid
cc_ctrl<-cc_ctrl[,-1]


##subset counts by vehicle 
CMC_AHR<-filter(meta1,vehicle == 'CMC_.5_%')
CORN_CARPXR<-filter(meta1,vehicle == 'CORN_OIL_100_%')
SALINE_Cytotoxic<-filter(meta1,vehicle == 'SALINE_100_%')

##Get names of filtered rows and filter merged counts by relevant vehicle column names
CMC_filter<-rownames(CMC_AHR)
CORN_filter<-rownames(CORN_CARPXR)
SALINE_Filter<-rownames(SALINE_Cytotoxic)

cnts_CMC<-select(cc_ctrl,"SRR1177998","SRR1178001","SRR1178003","SRR1178030","SRR1178040","SRR1178056")
cnts_CORN<-select(cc_ctrl, "SRR1177993", "SRR1177994" ,"SRR1177995", "SRR1178024", "SRR1178035", "SRR1178045")
cnts_SALINE<-select(cc_ctrl, "SRR1177966","SRR1177969","SRR1177970","SRR1178004","SRR1178006","SRR1178013")

#Reorder the rows of the meta data to match columns of counts 
indx <- match(colnames(cc_ctrl), rownames(meta1))
indx2<-indx[!is.na(indx)]
meta_x <- meta1[indx2,]
rownames(meta_x)



##rerun subset counts by vehicle post reorder
CMC_AHR<-filter(meta_x,vehicle == 'CMC_.5_%')
CORN_CARPXR<-filter(meta_x,vehicle == 'CORN_OIL_100_%')
SALINE_Cytotoxic<-filter(meta_x,vehicle == 'SALINE_100_%')


#check to confirm column names of counts file match row names for meta data
cc_mx<-select(cc_ctrl,rownames(meta_x))
all(rownames(meta_x) == colnames(cc_mx))
as.matrix(meta_x)
as.matrix(cc_mx)

# 4.3 create the DESeq object
dds_CMC_AHR <- DESeqDataSetFromMatrix(
  countData = cnts_CMC,
  colData = CMC_AHR,
  design= ~mode_of_action
)

dds_CORN_CARPXR <- DESeqDataSetFromMatrix(
  countData = cnts_CORN,
  colData = CORN_CARPXR,
  design= ~mode_of_action
)

dds_SALINE_Cytotoxic <- DESeqDataSetFromMatrix(
  countData = cnts_SALINE,
  colData = SALINE_Cytotoxic,
  design= ~mode_of_action
)

# relevel mode_of_action as factor - 
dds_CMC_AHR$mode_of_action <- relevel(dds_CMC_AHR$mode_of_action, ref='Control')
dds_CORN_CARPXR$mode_of_action <- relevel(dds_CORN_CARPXR$mode_of_action, ref='Control')
dds_SALINE_Cytotoxic$mode_of_action <- relevel(dds_SALINE_Cytotoxic$mode_of_action, ref='Control')

# run DESeq
dds_CMC_AHR <- DESeq(dds_CMC_AHR)
resCMC_AHR <- results(dds_CMC_AHR, contrast= c('mode_of_action','AhR','Control'))
resCMC_AHR <- lfcShrink(dds_CMC_AHR, coef=2)

dds_CORN_CARPXR <- DESeq(dds_CORN_CARPXR)
resCORN_CARPXR <- results(dds_CORN_CARPXR, contrast= c('mode_of_action','CAR/PXR','Control'))
resCORN_CARPXR <- lfcShrink(dds_CORN_CARPXR, coef=2)


dds_SALINE_Cytotoxic <- DESeq(dds_SALINE_Cytotoxic)
resSALINE_Cytotoxic <- results(dds_SALINE_Cytotoxic, contrast= c('mode_of_action','Cytotoxic','Control'))
resSALINE_Cytotoxic <- lfcShrink(dds_SALINE_Cytotoxic, coef=2)

##4.4
#Sort by adj p-value 
resCMC_AHR1 <- resCMC_AHR[order(resCMC_AHR$padj),]
resCORN_CARPXR1 <- resCORN_CARPXR[order(resCORN_CARPXR$padj),]
resSALINE_Cytotoxic1 <- resSALINE_Cytotoxic[order(resSALINE_Cytotoxic$padj),]

# write out DE results sorted by adjusted p-value
write.csv(resCMC_AHR1,'resCMC_AHR_deseq_results.csv')
write.csv(resCORN_CARPXR1,'resCORN_CARPXR_deseq_results.csv')
write.csv(resSALINE_Cytotoxic1,'resSALINE_Cytotoxic_deseq_results.csv')


#Report number of significant p-adjust <0.05
CMC_AHR_tab<-table(resCMC_AHR1$padj<0.05) 
CORN_CARPXR_tab<-table(resCORN_CARPXR1$padj<0.05)
SALINE_Cytotoxic_tab<-table(resSALINE_Cytotoxic1$padj<0.05)

all_sig<-bind_rows(CMC_AHR_tab,CORN_CARPXR_tab,SALINE_Cytotoxic_tab)
as.data.frame(all_sig)
crow = c("AHR & CMC .5%","CAR/PXR & CORN OIL 100%","Cytotoxic & SALINE 100%")
row.names(all_sig)<-crow


#4.5 Report top 10 genes by p-value
top_10_CMC_AHR<-read.csv('resCMC_AHR_deseq_results.csv')[1:10,]
top_10_CORN_CARPXR<-read.csv('resCORN_CARPXR_deseq_results.csv')[1:10,]
top_10_SALINE_Cytotoxic<-read.csv('resSALINE_Cytotoxic_deseq_results.csv')[1:10,]

#renaming X to Geneid
colnames(top_10_CMC_AHR)[colnames(top_10_CMC_AHR) == "X"] <- "Gene Id"
colnames(top_10_CORN_CARPXR)[colnames(top_10_CORN_CARPXR) == "X"] <- "Gene Id"
colnames(top_10_SALINE_Cytotoxic)[colnames(top_10_SALINE_Cytotoxic) == "X"] <- "Gene Id"



# 4.7 write out matrix of normalized counts
write.csv(counts(dds_CMC_AHR,normalized=TRUE),'CMC_AHR_deseq_norm_counts.csv')
write.csv(counts(dds_CORN_CARPXR,normalized=TRUE),'CORN__CARPXR_deseq_norm_counts.csv')
write.csv(counts(dds_SALINE_Cytotoxic,normalized=TRUE),'SALINE_Cytotoxic_deseq_norm_counts.csv')

##plot of foldchange 
hs_AHR_CMC<-hist(top_10_CMC_AHR$log2FoldChange, freq = TRUE, main = "Top 10 AHR Fold Change", xlab = "Fold Change (log2)", col ="pink", breaks = 15)

hs_CARPXR_CORN<-hist(top_10_CORN_CARPXR$log2FoldChange, col = "orange",freq = TRUE, main = "Top 10 CAR_PXR Fold Change", xlab = "Fold Change (log2)", breaks =15)

hs_Cytotoxic_SALINE<-hist(top_10_SALINE_Cytotoxic$log2FoldChange, col = "green",freq = TRUE, main = "Top 10 Cytotoxic Fold Change", xlab = "Fold Change (log2)",breaks = 15)

#scatter plots of p-value v foldchange

plot_AHR_CMC<-plot(x= top_10_CMC_AHR$log2FoldChange, y=  top_10_CMC_AHR$padj, main = "AHR Fold Change v p-value", xlab = "Fold Change (log2)", ylab="p-value (adj)", type ="p", col ="red")

plot_CARPXR_CORN<-plot(top_10_CORN_CARPXR$log2FoldChange, y= top_10_CORN_CARPXR$padj, main = "CAR_PXR Fold Change v p-value", xlab = "Fold Change (log2)", ylab="p-value(adj)", col = "blue")

plot_Cytotoxic_SALINE<-plot(top_10_SALINE_Cytotoxic$log2FoldChange, y= top_10_SALINE_Cytotoxic$padj, main = "Cytotoxic Change v p-value", xlab = "Fold Change (log2)", ylab="p-value(adj)", col = "green")

