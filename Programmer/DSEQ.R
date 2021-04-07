# 4.1 & 4.2  Setting up data for DESeq Analysis #
cnts <- as.matrix(cc_all2)

# sample information
info <- read.csv('allcounts.csv')
info2 = select(info, 1, 7:15)
info2[is.na(info2)] = 0
class(info2)
info2<-as.matrix(info2)
class(info2)

all(rownames(cnts) == colnames(meta1))

# 4.3 create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = cnts,
  colData = info2,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')

# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','HMGCOA','Control'))
res <- lfcShrink(dds, coef=2)

# write out DE results
write.csv(res,'example_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds,normalized=TRUE),'example_deseq_norm_counts.csv')