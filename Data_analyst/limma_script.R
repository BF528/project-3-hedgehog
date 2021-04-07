setwd("/projectnb/bf528/users/hedgehog/project_3/analyst")
library(limma)

# Reading in data ---------------------------------------------------------

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_2_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)


# BETA-NAPHTHOFLAVONE -----------------------------------------------------
## == CMC
# subset the full expression matrix to just those in this comparison
samples_betaNap <- samples$array_id[which(samples$chemical %in% c('BETA-NAPHTHOFLAVONE', "Control"))]
rma.subset <- rma[paste0('X',samples_betaNap)]

# construct a design matrix modeling treatment vs control for use by limma
design <- stats::model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','BETA-NAPHTHOFLAVONE')
  )
)
colnames(design) <- c('Control','BETA-NAPHTHOFLAVONE')

# run limma
fit <- limma::lmFit(rma.subset, design)
fit <- limma::eBayes(fit)
t1 <- limma::topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')



# ECONAZOLE ---------------------------------------------------------------
# == corn oil
# subset the full expression matrix to just those in this comparison
samples_econ <- samples$array_id[which(samples$chemical%in% c('ECONAZOLE', "Control"))]
rma.subset <- rma[paste0('X',samples_econ)]

# construct a design matrix modeling treatment vs control for use by limma
design <- stats::model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','ECONAZOLE')
  )
)
colnames(design) <- c('Control','ECONAZOLE')

# run limma
fit <- limma::lmFit(rma.subset, design)
fit <- limma::eBayes(fit)
t2 <- limma::topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')



# THIOACETAMIDE -----------------------------------------------------------
# == saline
# subset the full expression matrix to just those in this comparison
samples_thio <- samples$array_id[which(samples$chemical %in% c('THIOACETAMIDE', "Control"))]
rma.subset <- rma[paste0('X',samples_thio)]

design <- stats::model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','THIOACETAMIDE')
  )
)
colnames(design) <- c('Control','THIOACETAMIDE')

# run limma
fit <- limma::lmFit(rma.subset, design)
fit <- limma::eBayes(fit)
t3 <- limma::topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# Writing to csv ----------------------------------------------------------

# write out the results to file
# write.csv(t1,'limma_results_betaNap.csv')
# write.csv(t2,'limma_results_econ.csv')
# write.csv(t3,'limma_results_thio.csv')


# Getting DEs ---------------------------------------------------------
# Filtering for DEs
library(dplyr)
beta_genes <- t1 %>% 
  dplyr::filter(adj.P.Val<0.05)
econ_genes <- t2 %>% 
  dplyr::filter(adj.P.Val<0.05)
thio_genes <- t3 %>% 
  dplyr::filter(adj.P.Val<0.05)


# Deliverables -----------------------------------------------------------
## 1. "A report of the number of DE genes at p-adjust < 0.05"
beta_n <- nrow(beta_genes) # 3090 "genes"
econ_n <- nrow(econ_genes) # 5837 genes
thio_n <- nrow(thio_genes) # 8828 genes

beta_n0 <- nrow(t1) - beta_n
econ_n0 <- nrow(t2) - econ_n
thio_n0 <- nrow(t3) - thio_n

# report for number of DEGs
sig_g <- matrix(c(beta_n0, econ_n0, thio_n0, beta_n, econ_n, thio_n),ncol=2, nrow=3)
colnames(sig_g) <- c("FALSE", "TRUE")
rownames(sig_g) <- c('NAP', 'ECO', 'THI')

write.csv(sig_g, 'DEG_limma_report.csv')

## 2. "A table of the top 10 DE genes from each of your analyses"
beta_t10 <- beta_genes %>% 
  arrange(adj.P.Val) %>% 
  top_n(10)
econ_t10 <- econ_genes %>% 
  arrange(adj.P.Val) %>% 
  top_n(10)
thio_t10 <- thio_genes %>% 
  arrange(adj.P.Val) %>% 
  top_n(10)

# Writing them out
write.csv(beta_t10,'beta_t10.csv')
write.csv(econ_t10,'econ_t10.csv')
write.csv(thio_t10,'thio_t10.csv')

## 3. "Histograms of fold change values and scatter plots of fold change vs nominal-pvalue from the significant DE genes"
jpeg("histograms.jpeg")
par(mfrow=c(2,2), oma=c(0,0,2,0))
hist(beta_t10$logFC,
     main="Beta-Naphthoflavone",
     xlab="log(Fold Change)")
hist(econ_t10$logFC,
     main="Econazole",
     xlab="log(Fold Change)")
hist(thio_t10$logFC,
     main="Thioacetamide",
     xlab="log(Fold Change)")
mtext("Distribution of Fold Change Values", line = 0, outer=TRUE, cex=1.5)
dev.off()

jpeg("scatter_plots.jpeg")
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(beta_t10$logFC, beta_t10$P.Value,
     xlab = "log(FC)", ylab = "P-value",
     main = "Beta-Naphthoflavone")
plot(econ_t10$logFC, econ_t10$P.Value,
     xlab = "log(FC)", ylab = "P-value",
     main = "Econazole")
plot(thio_t10$logFC, thio_t10$P.Value,
     xlab = "log(FC)", ylab = "P-value",
     main = "Thioacetamide")
mtext("Scatter Plots of Fold Change vs P-value", line = 0, outer=TRUE, cex=1.5)
dev.off()