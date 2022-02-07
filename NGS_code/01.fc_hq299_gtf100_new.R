library(tximport)
library(edgeR)
library(mixOmics)
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
sample.pheno <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/299.IBD.pheno", what=character(0))
pca.pheno <- scan("/scratch/hpc46a05/gy_RNA/IBD/01.Samplelist/299.IBD.pheno", what=character(0))
names(sample.pheno) <- sample.ID

#condition - record - Pheno
sample.group <- factor(sample.pheno)
sample.group <- factor(substr(sample.pheno,1,5))
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)


#edgeR
cnt <- fc299
sample.DGElist <- DGEList(cnt, group=sample.group, genes=cnt[,0,drop=FALSE])
sample.filterExpr <- sample.DGElist[filterByExpr(sample.DGElist),]
sample.calNorm <- calcNormFactors(sample.filterExpr)
sample.estimateDisp <- estimateDisp(sample.calNorm,sample.design)
y <- sample.estimateDisp

#pca
sample.cpm <- cpm(sample.estimateDisp$counts,log=TRUE)
write.table(sample.cpm, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.IBD_gtf100_ensembl_cpm.csv", row.names=TRUE, quote=F, sep=",")
sample.cpm.trans <- t(sample.cpm)
sample.pca <- pca(sample.cpm.trans, ncomp=10, center=TRUE, scale=FALSE)

png(filename="/scratch/hpc46a05/gy_RNA/IBD/03.Plot/299.IBD.PCA.png")
plotIndiv(sample.pca, ncomp=c(1,2), ind.names=TRUE, group=pca.pheno, legend=TRUE, guide = "none")
dev.off()

#DEG
#[1] "CD_in" "CD_no" "UC_in" "UC_no"
DEG_fit <- glmQLFit(y, sample.design)

CDvsUC_in_DEG <- glmQLFTest(DEG_fit, contrast = c(1,0,-1,0))
CDvsUC_no_DEG <- glmQLFTest(DEG_fit, contrast = c(0,1,0,-1))
INvsNO_CD_DEG <- glmQLFTest(DEG_fit, contrast = c(1,-1,0,0))
INvsNO_UC_DEG <- glmQLFTest(DEG_fit, contrast = c(0,0,1,-1))
INvsNO_all_DEG <- glmQLFTest(DEG_fit, contrast = c(1,-1,1,-1))

CDvsUC_in <- decideTestsDGE(CDvsUC_in_DEG)
CDvsUC_no <- decideTestsDGE(CDvsUC_no_DEG)
INvsNO_CD <- decideTestsDGE(INvsNO_CD_DEG)
INvsNO_UC <- decideTestsDGE(INvsNO_UC_DEG)
INvsNO_all <- decideTestsDGE(INvsNO_all_DEG)

##Calc-FDR
CDvsUC_in_FDR <- p.adjust(CDvsUC_in_DEG$table$PValue, method="BH")
CDvsUC_no_FDR <- p.adjust(CDvsUC_no_DEG$table$PValue, method="BH")
INvsNO_CD_FDR <- p.adjust(INvsNO_CD_DEG$table$PValue, method="BH")
INvsNO_UC_FDR <- p.adjust(INvsNO_UC_DEG$table$PValue, method="BH")
INvsNO_all_FDR <- p.adjust(INvsNO_all_DEG$table$PValue, method="BH")

CDvsUC_in_result <- cbind(CDvsUC_in_DEG$table, CDvsUC_in_FDR)
CDvsUC_no_result <- cbind(CDvsUC_no_DEG$table, CDvsUC_no_FDR)
INvsNO_CD_result <- cbind(INvsNO_CD_DEG$table, INvsNO_CD_FDR)
INvsNO_UC_result <- cbind(INvsNO_UC_DEG$table, INvsNO_UC_FDR)
INvsNO_all_result <- cbind(INvsNO_all_DEG$table, INvsNO_all_FDR)

#DEG
write.table(y$counts, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.IBD_cnt.csv", row.names=TRUE, quote=F, sep=",")
write.table(CDvsUC_in_result, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.IBD_gtf100_CDvsUC_in_ensembl.csv", row.names=TRUE, quote=F, sep=",")
write.table(CDvsUC_no_result, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.IBD_gtf100_CDvsUC_no_ensembl.csv", row.names=TRUE, quote=F, sep=",")
write.table(INvsNO_CD_result, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.IBD_gtf100_INvsNO_CD_ensembl.csv", row.names=TRUE, quote=F, sep=",")
write.table(INvsNO_UC_result, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.IBD_gtf100_INvsNO_UC_ensembl.csv", row.names=TRUE, quote=F, sep=",")
write.table(INvsNO_all_result, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299.IBD_gtf100_INvsNO_all_ensembl.csv", row.names=TRUE, quote=F, sep=",")
