library(ggplot2)
library(RColorBrewer)

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
bx_cnts$out[8927]


##
hist()
