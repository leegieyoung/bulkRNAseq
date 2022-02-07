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
fc299 <- read.table("/scratch/hpc46a05/gy_RNA/IBD/02.featureCounts_result/Each_gene_id_P_M_countReadPairs_gtf100/merge/299.IBD.fc_P_M.merge", sep=' ', head=T)
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

#Find versus
CDvsUC_in <- c("CD_in", "UC_in")
INvsNO <- c("CD_in", "UC_in", "CD_no", "UC_no")
INvsNO_CD <- c("CD_in", "CD_no")
INvsNO_UC <- c("UC_in", "UC_no")

wCDvsUC_in <- which(sample.group %in% CDvsUC_in)
wINvsNO <- which(sample.group %in% INvsNO)
wINvsNO_CD <- which(sample.group %in% INvsNO_CD)
wINvsNO_UC <- which(sample.group %in% INvsNO_UC)

#anno_function
gtf <- read.table("/scratch/hpc46a05/gy_RNA/Reference/gtf/anno_Homo_sapiens.GRCh38.100.gtf", sep=" ", head=T)

cnt <- fc299

#logCPM

#edgeR
cnt.lcpm <- cpm(cnt, log=TRUE)
y <- DGEList(cnt, group=sample.group, genes=rownames(cnt.lcpm),)

#Ensembl add function
idfound <- rownames(y) %in% mappedRkeys(org.Hs.egENSEMBL)
y <- y[idfound,]
m <- match(y$genes$genes, gtf$ensembl_id)
Func <- gtf$function.[m]
y$genes <- cbind(y$genes, Func)

#Ensembl to Entrezid
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

#Entrezid -> Gene Symbol
egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(rownames(y_uniq), egSYMBOL$gene_id)
y_uniq$genes$genes <- egSYMBOL$symbol[m]

keep <- filterByExpr(y_uniq)
y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
#y_uniq <- y_uniq[filterByExpr(y_uniq),]
y_uniq <- calcNormFactors(y_uniq)
y_uniq <- estimateDisp(y_uniq,sample.design)
y_uniq.lcpm <- cpm(y_uniq, log=TRUE)

#glmQLF = F-test
Ffit <- glmQLFit(y_uniq, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=c(1,-1,0,0))
png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/INvsNO_CD/01.Plot/glmQLF.png", width=60, height=50, units="cm", res=200)
plotMD(Ftest)
abline(h=c(-1, 1), col="blue")
dev.off()

print(Ftest$comparison)

#decideTests
dt <- (decideTests(Ftest, p.value=0.05, adjust.method="fdr", lfc=1))
w <- which(dt==0)
qcFtest <- Ftest[-w,]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/INvsNO_CD/01.Plot/QC_glmQLF.png", width=60, height=50, units="cm", res=200)
plotMD(qcFtest)
abline(h=c(-1, 1), col="blue")
dev.off()

#dt Gene Heatmap, This purpose is that find cpm of not value.
w <- which(rownames(y_uniq) %in% rownames(qcFtest$genes))
y_uniq_dt <- y_uniq[w,]
y_uniq_dt.lcpm <- cpm(y_uniq_dt, log=T)
y_uniq_dt.lcpm <- y_uniq_dt.lcpm[,wINvsNO_CD]
heatmap.pheno <- sample.pheno[wINvsNO_CD]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/INvsNO_CD/01.Plot/uniq_dt_heatmap.png", width=60, height=150, units="cm", res=200)
col <- list(IBD = c("CD_in" = "red", "CD_no"="green"), INvsNO = c("CD_in"="purple", "CD_no"="orange"))
anno <- HeatmapAnnotation(IBD = heatmap.pheno ,INvsNO=heatmap.pheno, col=col)
Heatmap(y_uniq_dt.lcpm, top_annotation = anno, column_names_gp=grid::gpar(fontsize=12), row_names_gp = grid::gpar(fontsize = 12), width=unit(45, "cm"), height=unit(135,"cm"))
dev.off()

#Check contrast values.
print("Check contrast values.")
colnames(Ffit)
#[1] "CD_in" "CD_no" "UC_in" "UC_no"
table <- cbind(qcFtest$table, qcFtest$genes)
write.table(table, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/INvsNO_CD/02.Table/genes.csv", quote=F, col.names=T, row.names=T,sep=",")

#pathway-go-kegg
go.CDvsUC_in <- goana(qcFtest, species="Hs")
keg.CDvsUC_in <- kegga(qcFtest, species="Hs")
cutoff.go <- subset(go.CDvsUC_in, go.CDvsUC_in$P.Up < 0.05 | go.CDvsUC_in$P.Down < 0.05)
cutoff.go$Term <- gsub(" ", "_", cutoff.go$Term)
cutoff.go$Term <- gsub(",", "|", cutoff.go$Term)
go.up.down <- cbind(rownames(cutoff.go), cutoff.go$Term, cutoff.go$P.Up, cutoff.go$P.Down)
write.table(go.up.down,"/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/INvsNO_CD/03.Pathway/go.up.down.csv",quote=F, col.names=T, row.names=T,sep=",")

cutoff.keg <- subset(keg.CDvsUC_in, keg.CDvsUC_in$P.Up < 0.05 | keg.CDvsUC_in$P.Down < 0.05)
cutoff.keg$Pathway <- gsub(",", "|", cutoff.keg$Pathway)
cutoff.keg$Pathway <- gsub(" ", "_", cutoff.keg$Pathway)
keg.up.down <- cbind(rownames(cutoff.keg), cutoff.keg$Pathway, cutoff.keg$P.Up, cutoff.keg$P.Down)
write.table(keg.up.down,"/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/INvsNO_CD/03.Pathway/kegg.up.down.csv",quote=F, col.names=T, row.names=T,sep=",")


#gene
#tr <- glmTreat(fit, contrast=c(1,0,-1,0), lfc=log2(2))
#tr_b <- cbind(tr$table, tr$genes)
#tr_cutoff <- subset(tr_b, tr_b$PValue < 0.05)

#CDvsUC_in <- decideTestsDGE(CDvsUC_in_DEG)
