library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(repr)
library(statmod)
library(GO.db)
library(ComplexHeatmap)
library(org.Hs.eg.db)

#=============
#Scan
sample.ID <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/454.IBD.sample", what=character(0))
sample.raw <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/454.IBD.pheno", what=character(0))
pca.pheno <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/454.IBD.pheno", what=character(0))
names(sample.raw) <- sample.ID
sample.pheno <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/454.IBD.pheno", what=character(0))
sample.pheno <- substr(sample.pheno,1,5)

#condition - record - Pheno
sample.group <- factor(sample.raw)
sample.group <- factor(substr(sample.raw,1,5))
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)

#299 Totalseq
Totalseq.ID <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/299.IBD.sample", what=character(0))
w <- which(sample.ID %in% Totalseq.ID)
Totalseq.raw <- sample.raw[w]
Totalseq.group <- factor(Totalseq.raw)
Totalseq.group <- factor(substr(Totalseq.raw,1,5))
Totalseq.design <- model.matrix(~0 + Totalseq.group)
colnames(Totalseq.design) <- levels(Totalseq.group)

#155 Truseq is not Totalseq in 454samples
Truseq.ID <- sample.ID[-w]
Truseq.raw <- sample.raw[-w]
Truseq.group <- factor(Truseq.raw)
Truseq.group <- factor(substr(Truseq.raw,1,5))
Truseq.design <- model.matrix(~0 + Truseq.group)
colnames(Truseq.design) <- levels(Truseq.group)


cx299 <- read.table("/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/02.Cibersortx/299/CIBERSORTx_299_merge.csv", sep=",", head=T)
rownames(cx299) <- cx299[,1]
cx299 <- cx299[,-1]
cx299 <- cx299[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)]

plotIndiv(cx299.splsda, comp=1:2, group=Totalseq.group, ellipse = TRUE)


