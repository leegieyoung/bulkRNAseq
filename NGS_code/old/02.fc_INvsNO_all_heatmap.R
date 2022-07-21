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

#Do change select_miRNA
select_miRNA <- scan("/scratch/hpc46a05/gy_RNA/IBD/04.CSV/test/select_miRNA", what=character(0))

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


#edgeR
cnt <- fc299
sample.DGElist <- DGEList(cnt, group=sample.group, genes=cnt[,0,drop=FALSE])
sample.filterExpr <- sample.DGElist[filterByExpr(sample.DGElist),]
sample.calNorm <- calcNormFactors(sample.filterExpr)
sample.estimateDisp <- estimateDisp(sample.calNorm,sample.design)
y <- sample.estimateDisp

#logCPM
y.lcpm <- cpm(y,log=TRUE)

INvsNO_fit <- glmQLFit(y, sample.design)
colnames(INvsNO_fit)

#[1] "CD_in" "CD_no" "UC_in" "UC_no"

#heatmap
#IBD and INvsNO
INvsNO_DEG <- glmQLFTest(INvsNO_fit, contrast = c(1,-1,1,-1))
is.de <- decideTestsDGE(INvsNO_DEG)
INvsNO_DEGlist <- is.de[!(is.de == "0" ), ]
INvsNO_lcpm <- y.lcpm[rownames(INvsNO_DEGlist),]

highly_variable_lcpm <- INvsNO_lcpm[select_miRNA,]
col <- list(IBD = c("CD_in" = "red", "UC_no"="blue", "CD_no"="green", "UC_in"="yellow"), INvsNO = c("CD_in"="purple", "UC_in"="purple", "CD_no"="orange", "UC_no"="orange"))
anno <- HeatmapAnnotation(IBD = sample.pheno ,INvsNO=sample.pheno, col=col)

png(filename="/scratch/hpc46a05/gy_RNA/IBD/03.Plot/299.heatmap.png", width=48, height=20, units="cm", res=1000)
Heatmap(highly_variable_lcpm, top_annotation = anno, column_names_gp=grid::gpar(fontsize=6), row_names_gp = grid::gpar(fontsize = 6), width=unit(36, "cm"), height=unit(16,"cm"))
dev.off()

