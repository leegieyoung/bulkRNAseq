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

#hisat-featurecounts result load
fc299 <- read.table("/scratch/hpc46a05/gy_RNA/IBD/02.featureCounts_result/Each_gene_id_P_M_countReadPairs_gtf100/merge/299.IBD.fc_P_M.merge", sep=' ',head=T)
fc299.ensembl_id <- fc299[,1]
fc299 <- fc299[,-1]

#rownames
rownames(fc299) <- fc299.ensembl_id

#Simplify smaples name & pheno : dot -> hyphen
colnames(fc299) <- gsub(".TR.sorted.bam","",colnames(fc299))
colnames(fc299) <- gsub("[.]","-",colnames(fc299))

#Scan
sample.ID <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/299.IBD.sample", what=character(0))
sample.raw <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/299.IBD.pheno", what=character(0))
pca.pheno <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/299.IBD.pheno", what=character(0))
names(sample.raw) <- sample.ID
sample.pheno <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/299.IBD.pheno", what=character(0))
sample.pheno <- substr(sample.pheno,1,5)

#condition - record - Pheno
sample.group <- factor(sample.raw)
sample.group <- factor(substr(sample.raw,1,5))
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)

cnt <- fc299

#logCPM
cnt.lcpm <- cpm(cnt, log=TRUE)
write.table(cnt.lcpm, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.IBD_gtf100_INvsNO_all_ensembl/cpm.txt", quote=F, sep="\t")

#edgeR
y <- DGEList(cnt, group=sample.group, genes=rownames(cnt.lcpm),)

#Ensembl to Entrezid
idfound <- rownames(y) %in% mappedRkeys(org.Hs.egENSEMBL)
y <- y[idfound,]

egENTREZID <- toTable(org.Hs.egENSEMBL)
m <- match(rownames(y), egENTREZID$ensembl_id)
y$genes$genes <- egENTREZID$gene_id[m]

#Uniq Entrezid
o <- order(rowSums(y$counts), decreasing=TRUE)
y_uniq <- y[o,]
d <- duplicated(y_uniq$genes$genes)
y_uniq <- y_uniq[!d,]
rownames(y_uniq$counts) <- rownames(y_uniq$genes) <- y_uniq$genes$genes
y_uniq$genes$genes <- NULL

cnt <- fc299
keep <- filterByExpr(y_uniq)
y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
#y_uniq <- y_uniq[filterByExpr(y_uniq),]
y_uniq <- calcNormFactors(y_uniq)
y_uniq <- estimateDisp(y_uniq,sample.design)

INvsNO_fit <- glmQLFit(y_uniq, sample.design)
colnames(INvsNO_fit)
#[1] "CD_in" "CD_no" "UC_in" "UC_no"

#CD_Match cmap
cmap <- scan("/scratch/hpc46a05/gy_RNA/Reference/cMAP/entrezid.txt",what=character(0))

print("===============================")
print("Do reading Cmap annotation file")
print("===============================")

#CD
DEG <- glmQLFTest(INvsNO_fit, contrast = c(1,-1,0,0))

print("===================")
print("Complete glmQLFTest")
print("===================")

#CD_Up
DEG_U <- DEG$table[DEG$table$logFC > 1 & DEG$table$PValue < 0.05,]
m <- match(cmap,rownames(DEG_U[order(-DEG_U$logFC),]))
DEG_U <- na.omit(DEG_U[m,])
cmap_up <- head(DEG_U[order(-DEG_U$logFC),],50)
write.table(cmap_up, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.cmap/CD_up.csv", sep=",", col.names=T, row.names=T, quote=F)

#CD_Down
DEG_D <- DEG$table[DEG$table$logFC < -1 & DEG$table$PValue < 0.05,]
m <- match(cmap,rownames(DEG_D[order(DEG_D$logFC),]))
DEG_D <- na.omit(DEG_D[m,])
cmap_dn <- head(DEG_D[order(DEG_D$logFC),],50)
write.table(cmap_up, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.cmap/CD_down.csv", sep=",", col.names=T, row.names=T, quote=F)

#UC
DEG <- glmQLFTest(INvsNO_fit, contrast = c(0,0,1,-1))
print("===================")
print("Complete glmQLFTest")
print("===================")

#UC_Up
DEG_U <- DEG$table[DEG$table$logFC > 1 & DEG$table$PValue < 0.05,]
m <- match(cmap,rownames(DEG_U[order(-DEG_U$logFC),]))
DEG_U <- na.omit(DEG_U[m,])
cmap_up <- head(DEG_U[order(-DEG_U$logFC),],50)
write.table(cmap_up, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.cmap/UC_up.csv", sep=",", col.names=T, row.names=T, quote=F)

#UC_Down
DEG_D <- DEG$table[DEG$table$logFC < -1 & DEG$table$PValue < 0.05,]
m <- match(cmap,rownames(DEG_D[order(DEG_D$logFC),]))
DEG_D <- na.omit(DEG_D[m,])
cmap_dn <- head(DEG_D[order(DEG_D$logFC),],50)
write.table(cmap_up, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.cmap/UC_down.csv", sep=",", col.names=T, row.names=T, quote=F)

