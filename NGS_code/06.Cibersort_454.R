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
fc454 <- read.table("/scratch/hpc46a05/gy_RNA/IBD/02.featureCounts_result/Each_gene_id_P_M_countReadPairs_gtf100/merge/454.IBD.fc_P_M.merge", sep=' ',head=T)
fc454.ensembl_id <- fc454[,1]
fc454 <- fc454[,-1]

#rownames
rownames(fc454) <- fc454.ensembl_id

#Simplify smaples name & pheno : dot -> hyphen
colnames(fc454) <- gsub(".TR.sorted.bam","",colnames(fc454))
colnames(fc454) <- gsub("[.]","-",colnames(fc454))

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
fc299 <- fc454[w]
Totalseq.raw <- sample.raw[w]
Totalseq.group <- factor(Totalseq.raw)
Totalseq.group <- factor(substr(Totalseq.raw,1,5))
Totalseq.design <- model.matrix(~0 + Totalseq.group)
colnames(Totalseq.design) <- levels(Totalseq.group)

#155 Truseq is not Totalseq in 454samples
fc155 <- fc454[-w]
Truseq.ID <- sample.ID[-w]
Truseq.raw <- sample.raw[-w]
Truseq.group <- factor(Truseq.raw)
Truseq.group <- factor(substr(Truseq.raw,1,5))
Truseq.design <- model.matrix(~0 + Truseq.group)
colnames(Truseq.design) <- levels(Truseq.group)



#===============================================================
#fc454 and DGElist and annotation(Gene Symbol)
cnt <- fc454
cnt.lcpm <- cpm(cnt, log=TRUE)

#edgeR
y <- DGEList(cnt, group=sample.group, genes=rownames(cnt.lcpm),)
#Gene : 60683
#Ensembl to Entrezid
idfound <- rownames(y) %in% mappedRkeys(org.Hs.egENSEMBL)
y <- y[idfound,]
egSYMBOL <- toTable(org.Hs.egSYMBOL)
#dim : 63901
egENTREZID <- toTable(org.Hs.egENSEMBL)
#dim : 38849
m <- match(rownames(y), egENTREZID$ensembl_id)
#첫번째 벡터의 인수가 두번째 벡터의 인수의 몇번째에 있는지를 알려줌
y$genes$genes <- egENTREZID$gene_id[m]
#Uniq Entrezid
o <- order(rowSums(y$counts), decreasing=TRUE)
y_uniq <- y[o,]
d <- duplicated(y_uniq$genes$genes)
y_uniq <- y_uniq[!d,]
rownames(y_uniq$counts) <- rownames(y_uniq$genes) <- y_uniq$genes$genes

#find Symbol
m <- match(rownames(y_uniq), egSYMBOL$gene_id)
rownames(y_uniq) <- egSYMBOL$symbol[m]
o <- order(rowSums(y_uniq$counts), decreasing=TRUE)
y_uniq <- y_uniq[o,]
d <- duplicated(y_uniq$genes$genes)
y_uniq <- y_uniq[!d,]

totalseq_uniq <- y_uniq[,w]
truseq_uniq <- y_uniq[,-w]

#filterByExpr and calcNormFactors
#totalseq
keep <- filterByExpr(totalseq_uniq)
totalseq_uniq <- totalseq_uniq[keep, , keep.lib.sizes=FALSE]
totalseq_uniq <- calcNormFactors(totalseq_uniq)

#truseq
keep <- filterByExpr(truseq_uniq)
truseq_uniq <- truseq_uniq[keep, , keep.lib.sizes=FALSE]
truseq_uniq <- calcNormFactors(truseq_uniq)

#Genes of intersection and extract intersection genes
m <- match(rownames(truseq_uniq), rownames(totalseq_uniq))
m <- na.omit(m)
totalseq_uniq <- (totalseq_uniq[m,])

m <- match(rownames(totalseq_uniq), rownames(truseq_uniq))
m <- na.omit(m)
truseq_uniq <- (truseq_uniq[m,])

totalseq_uniq.lcpm <- cpm(totalseq_uniq, log=TRUE)
truseq_uniq.lcpm <- cpm(truseq_uniq, log=TRUE)

/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/intersect/01.edgeR
#total 299
write.table(totalseq_uniq.lcpm, "/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/299.IBD.txt", sep="\t", row.names=T, col.names=T, quote=F)
write.table(truseq_uniq.lcpm, "/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/155.IBD.txt", sep="\t", row.names=T, col.names=T, quote=F)

#UC_in
w <- which(Totalseq.group == "UC_in")
n <- names(Totalseq.group[w])
UC_in <- totalseq_uniq.lcpm[,n]
write.table(UC_in, "/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/299.UC_in.txt", sep="\t", row.names=T, col.names=T, quote=F)

w <- which(Truseq.group == "UC_in")
n <- names(Truseq.group[w])
UC_in <- truseq_uniq.lcpm[,n]
write.table(UC_in, "/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/155.UC_in.txt", sep="\t", row.names=T, col.names=T, quote=F)

#UC_no
w <- which(Totalseq.group == "UC_no")
n <- names(Totalseq.group[w])
UC_no <- totalseq_uniq.lcpm[,n]
write.table(UC_no, "/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/299.UC_no.txt", sep="\t", row.names=T, col.names=T, quote=F)

w <- which(Truseq.group == "UC_no")
n <- names(Truseq.group[w])
UC_no <- truseq_uniq.lcpm[,n]
write.table(UC_no, "/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/155.UC_no.txt", sep="\t", row.names=T, col.names=T, quote=F)

#CD_in
w <- which(Totalseq.group == "CD_in")
n <- names(Totalseq.group[w])
CD_in <- totalseq_uniq.lcpm[,n]
write.table(CD_in, "/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/299.CD_in.txt", sep="\t", row.names=T, col.names=T, quote=F)

w <- which(Truseq.group == "CD_in")
n <- names(Truseq.group[w])
CD_in <- truseq_uniq.lcpm[,n]
write.table(CD_in, "/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/155.CD_in.txt", sep="\t", row.names=T, col.names=T, quote=F)

#CD_no
w <- which(Totalseq.group == "CD_no")
n <- names(Totalseq.group[w])
CD_no <- totalseq_uniq.lcpm[,n]
write.table(CD_no, "/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/299.CD_no.txt", sep="\t", row.names=T, col.names=T, quote=F)

w <- which(Truseq.group == "CD_no")
n <- names(Truseq.group[w])
CD_no <- truseq_uniq.lcpm[,n]
write.table(CD_no, "/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/intersect/01.edgeR/155.CD_no.txt", sep="\t", row.names=T, col.names=T, quote=F)

