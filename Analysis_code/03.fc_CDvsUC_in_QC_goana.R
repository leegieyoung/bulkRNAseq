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
library(mixOmics)

print("Loading raw data...")
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

#For contrast
#[1] "CD_in" "CD_no" "UC_in" "UC_no"
my.contrasts <- makeContrasts(cCDvsUC_in=CD_in-UC_in, cCDvsUC_no=CD_no-UC_no, cINvsNO_CD=CD_in-CD_no, cINvsNO_UC=UC_in-UC_no,cINvsNO_all=CD_in+UC_in-CD_no-UC_no, levels=sample.design)

#Find versus
CDvsUC_in <- c("CD_in", "UC_in")
CDvsUC_no <- c("CD_no", "UC_no")
INvsNO_all <- c("CD_in", "UC_in", "CD_no", "UC_no")
INvsNO_CD <- c("CD_in", "CD_no")
INvsNO_UC <- c("UC_in", "UC_no")

wCDvsUC_in <- which(sample.group %in% CDvsUC_in)
wCDvsUC_no <- which(sample.group %in% CDvsUC_no)
wINvsNO_all <- which(sample.group %in% INvsNO_all)
wINvsNO_CD <- which(sample.group %in% INvsNO_CD)
wINvsNO_UC <- which(sample.group %in% INvsNO_UC)

CDvsUC_in.raw <- sample.raw[wCDvsUC_in]
CDvsUC_in.group <- factor(CDvsUC_in.raw)
CDvsUC_in.group <- factor(substr(CDvsUC_in.raw,1,5))

CDvsUC_no.raw <- sample.raw[wCDvsUC_no]
CDvsUC_no.group <- factor(CDvsUC_no.raw)
CDvsUC_no.group <- factor(substr(CDvsUC_no.raw,1,5))

INvsNO_all.raw <- sample.raw[wINvsNO_all]
INvsNO_all.group <- factor(INvsNO_all.raw)
INvsNO_all.group <- factor(substr(INvsNO_all.raw,1,5))

INvsNO_CD.raw <- sample.raw[wINvsNO_CD]
INvsNO_CD.group <- factor(INvsNO_CD.raw)
INvsNO_CD.group <- factor(substr(INvsNO_CD.raw,1,5))

INvsNO_UC.raw <- sample.raw[wINvsNO_UC]
INvsNO_UC.group <- factor(INvsNO_UC.raw)
INvsNO_UC.group <- factor(substr(INvsNO_UC.raw,1,5))

#For heatmap
col_CDvsUC_in <- list(IBD = c("CD_in" = "red", "UC_in"="yellow"), INvsNO = c("CD_in"="purple", "UC_in"="purple"))
col_CDvsUC_no <- list(IBD = c("CD_no" = "red", "UC_no"="yellow"), INvsNO = c("CD_no"="purple", "UC_no"="purple"))
col_INvsNO_all <- list(IBD = c("CD_in" = "red", "UC_no"="blue", "CD_no"="green", "UC_in"="yellow"), INvsNO = c("CD_in"="purple", "UC_in"="purple", "CD_no"="orange", "UC_no"="orange"))
col_INvsNO_CD <- list(IBD = c("CD_in" = "red",  "CD_no"="green"), INvsNO = c("CD_in"="purple", "CD_no"="orange"))
col_INvsNO_UC <- list(IBD = c("UC_no"="blue", "UC_in"="yellow"), INvsNO = c("UC_in"="purple", "UC_no"="orange"))

#anno_function
gtf <- read.table("/scratch/hpc46a05/gy_RNA/Reference/gtf/anno_Homo_sapiens.GRCh38.100.gtf", sep=" ", head=T)

print("End of loading raw data")

cnt <- fc299

print("Do edgeR and Annotation using org.Hs")
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

#Remove LOC gene
remove <- grep("^LOC", y_uniq$genes$genes, value=T)
y_uniq <- y_uniq[which(!y_uniq$genes$genes %in% remove),]

keep <- filterByExpr(y_uniq)
y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
#y_uniq <- y_uniq[filterByExpr(y_uniq),]
y_uniq <- calcNormFactors(y_uniq)
y_uniq <- estimateDisp(y_uniq,sample.design)
y_uniq.lcpm <- cpm(y_uniq, log=TRUE)

#glmQLF = F-test
Ffit <- glmQLFit(y_uniq, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=my.contrasts[,"cCDvsUC_in"])
png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/01.Plot/glmQLF.png", width=60, height=50, units="cm", res=200)
plotMD(Ftest)
abline(h=c(-1, 1), col="blue")
dev.off()

print(Ftest$comparison)

#decideTests
dt <- (decideTests(Ftest, p.value=0.05, adjust.method="fdr", lfc=1))
w <- which(dt==0)
qcFtest <- Ftest[-w,]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/01.Plot/QC_glmQLF.png", width=60, height=50, units="cm", res=200)
plotMD(qcFtest)
abline(h=c(-1, 1), col="blue")
dev.off()

print("Draw Heatmap")
#dt Gene Heatmap, This purpose is that find cpm of not value.
w <- which(rownames(y_uniq) %in% rownames(qcFtest$genes))
y_uniq_dt <- y_uniq[w,]
y_uniq_dt.lcpm <- cpm(y_uniq_dt, log=T)
y_uniq_dt.lcpm <- y_uniq_dt.lcpm[,wCDvsUC_in]
heatmap.pheno <- sample.pheno[wCDvsUC_in]

#Scale center = T, scale = F
test <- y_uniq_dt.lcpm
for(i in 1:length(rownames(test))){
test[i,] <- scale(test[i,], center=TRUE, scale=FALSE)
}
png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/01.Plot/uniq_dt_centerT_heatmap.png", width=60, height=160, units="cm", res=200)
anno <- HeatmapAnnotation(
	IBD = heatmap.pheno ,INvsNO=heatmap.pheno, col=col_CDvsUC_in, 
	simple_anno_size = unit(2.5, "cm"),height = unit(2.5, "cm"), 
	annotation_name_rot = 45,
	annotation_name_gp=gpar(fontsize=50), 
	annotation_legend_param=list(
		IBD = list(title_gp=gpar(fontsize=40), labels_gp=gpar(fontsize=35)),
		INvsNO= list(title_gp=gpar(fontsize=40), labels_gp=gpar(fontsize=35))	
	)
)

Heatmap(test, top_annotation = anno, 
	column_names_gp=grid::gpar(fontsize=12), row_names_gp = grid::gpar(fontsize = 12), 
	width=unit(45, "cm"), height=unit(135,"cm"),
	heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 30))
)
dev.off()

png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/01.Plot/uniq_dt_heatmap.png", width=60, height=160, units="cm", res=200)
Heatmap(y_uniq_dt.lcpm, top_annotation = anno,                                                                                                                                                 column_names_gp=grid::gpar(fontsize=12), row_names_gp = grid::gpar(fontsize = 12),
        width=unit(45, "cm"), height=unit(135,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 30))
)
dev.off()

#Heatmap NLRP3, NOD2
NLRP3 <- which(y_uniq_dt$genes$genes %in% "NLRP3")
NOD2 <- which(y_uniq_dt$genes$genes %in% "NOD2")
sig_gene <- c(NLRP3, NOD2)
y_sig_gene_dt.lcpm <- y_uniq_dt.lcpm[sig_gene,]
rownames(y_sig_gene_dt.lcpm) <- c("NLRP3", "NOD2")
#Scale center = T, scale = F

test <- y_sig_gene_dt.lcpm
for(i in 1:length(rownames(test))){
test[i,] <- scale(test[i,], center=TRUE, scale=FALSE)
}
png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/01.Plot/sig_gene_dt_centerT_heatmap.png", width=120, height=40, units="cm", res=200)
anno <- HeatmapAnnotation(
        IBD = heatmap.pheno ,INvsNO=heatmap.pheno, col=col_CDvsUC_in,
        simple_anno_size = unit(2.5, "cm"),height = unit(2.5, "cm"),
        annotation_name_rot = 45,
        annotation_name_gp=gpar(fontsize=50),
        annotation_legend_param=list(
                IBD = list(title_gp=gpar(fontsize=40), labels_gp=gpar(fontsize=35)),
                INvsNO= list(title_gp=gpar(fontsize=40), labels_gp=gpar(fontsize=35))
        )
)

Heatmap(test, top_annotation = anno,
        column_names_gp=grid::gpar(fontsize=20), row_names_gp = grid::gpar(fontsize = 40),
        width=unit(100, "cm"), height=unit(10,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 30))
)
dev.off()

png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/01.Plot/sig_gene_dt_heatmap.png", width=120, height=40, units="cm", res=200)
Heatmap(y_sig_gene_dt.lcpm, top_annotation = anno,
	column_names_gp=grid::gpar(fontsize=20), row_names_gp = grid::gpar(fontsize = 40),
	width=unit(100, "cm"), height=unit(10,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 30))
)
dev.off()

print("Draw PCA")
#PCA
y_uniq_dt.lcpm.t <- t(y_uniq_dt.lcpm)
y_uniq_dt.pca <- pca(y_uniq_dt.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/01.Plot/uniq_dt_PCA.png")
plotIndiv(y_uniq_dt.pca, ncomp=c(1,2), ind.names=TRUE, group=heatmap.pheno, legend=TRUE, guide = "none", title="PCA")
dev.off()

print("Draw sPLS-DA")
#sPLS-DA
y_uniq_dt.lcpm.t <- t(y_uniq_dt.lcpm[,which(colnames(y_uniq_dt.lcpm) %in% sample.ID[wCDvsUC_in])])
y_uniq_dt.splsda <- splsda(y_uniq_dt.lcpm.t, CDvsUC_in.group, ncomp=10)
png(filename="/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/01.Plot/uniq_dt_splsda.png")
plotIndiv(y_uniq_dt.splsda, comp = 1:2, 
	group = CDvsUC_in.group, ind.names = FALSE,
	ellipse = TRUE, 
	legend = TRUE,
	title="sPLS-DA"
)
background = background.predict(y_uniq_dt.splsda, comp.predicted=2, dist = "max.dist")
dev.off()

print("Write genes")
#Check contrast values.
print("Check contrast values.")
colnames(Ffit)
#[1] "CD_in" "CD_no" "UC_in" "UC_no"
table <- cbind(qcFtest$table, qcFtest$genes)
GOnet <- cbind(qcFtest$genes$genes ,qcFtest$table$logFC)
David <- cbind(rownames(qcFtest$genes) ,qcFtest$table$logFC)
write.table(table, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/02.Table/genes.csv", quote=F, col.names=T, row.names=T,sep=",")
write.table(GOnet, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/02.Table/GOnet.txt", quote=F, col.names=F, row.names=F,sep="\t")
write.table(David, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/02.Table/David.txt", quote=F, col.names=F, row.names=F,sep="\t")
print("Draw pathway")
#pathway-go-kegg
go.CDvsUC_in <- goana(qcFtest, species="Hs")
keg.CDvsUC_in <- kegga(qcFtest, species="Hs")
cutoff.go <- subset(go.CDvsUC_in, go.CDvsUC_in$P.Up < 0.05 | go.CDvsUC_in$P.Down < 0.05)
cutoff.go$Term <- gsub(" ", "_", cutoff.go$Term)
cutoff.go$Term <- gsub(",", "|", cutoff.go$Term)
go.up.down <- cbind(rownames(cutoff.go), cutoff.go$Term, cutoff.go$P.Up, cutoff.go$P.Down)
write.table(go.up.down,"/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/03.Pathway/go.up.down.csv",quote=F, col.names=F, row.names=F, sep=",")

cutoff.keg <- subset(keg.CDvsUC_in, keg.CDvsUC_in$P.Up < 0.05 | keg.CDvsUC_in$P.Down < 0.05)
cutoff.keg$Pathway <- gsub(",", "|", cutoff.keg$Pathway)
cutoff.keg$Pathway <- gsub(" ", "_", cutoff.keg$Pathway)
keg.up.down <- cbind(rownames(cutoff.keg), cutoff.keg$Pathway, cutoff.keg$P.Up, cutoff.keg$P.Down)
write.table(keg.up.down,"/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/03.Pathway/kegg.up.down.csv",quote=F, col.names=F, row.names=F, sep=",")

#go_term
go_cg <- read.table("/scratch/hpc46a05/gy_RNA/Reference/GOterm/go_categories.txt", sep=",",head =F)
go_namespace <- c("go","namespace")
colnames(go_cg) <- go_namespace
w <- which(go_cg$namespace %in% "biological_process")
bp <- go_cg[w,]
w <- which(go_cg$namespace %in% "molecular_function")
mf <- go_cg[w,]
w <- which(go_cg$namespace %in% "cellular_component")
cc <- go_cg[w,]
go_p <- c("GO","P")

go.up <- cutoff.go[which(cutoff.go$P.Up < 0.05),]
go.up <- cbind(rownames(go.up), go.up$Term, go.up$P.Up)

go.down <- cutoff.go[which(cutoff.go$P.Down < 0.05),]
go.down <- cbind(rownames(go.down), go.down$Term, go.down$P.Down)

go.up.bp <- go.up[which(go.up[,1] %in% bp[,1]),]
write.table(go.up.bp, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/03.Pathway/go.up.bp.txt",quote=F, col.names=F, row.names=F, sep="\t")
go.up.cc <- go.up[which(go.up[,1] %in% cc[,1]),]
write.table(go.up.cc, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/03.Pathway/go.up.cc.txt",quote=F, col.names=F, row.names=F, sep="\t")
go.up.mf <- go.up[which(go.up[,1] %in% mf[,1]),]
write.table(go.up.mf, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/03.Pathway/go.up.mf.txt",quote=F, col.names=F, row.names=F, sep="\t")

go.down.bp <- go.down[which(go.down[,1] %in% bp[,1]),]
write.table(go.down.bp, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/03.Pathway/go.down.bp.txt",quote=F, col.names=F, row.names=F, sep="\t")
go.down.cc <- go.down[which(go.down[,1] %in% cc[,1]),]
write.table(go.down.cc, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/03.Pathway/go.down.cc.txt",quote=F, col.names=F, row.names=F, sep="\t")
go.down.mf <- go.down[which(go.down[,1] %in% mf[,1]),]
write.table(go.down.mf, "/scratch/hpc46a05/gy_RNA/IBD/04.CSV/299/CDvsUC_in/03.Pathway/go.down.mf.txt",quote=F, col.names=F, row.names=F, sep="\t")
