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

#====================================
#299
#CD_in
CD_in <- read.table("/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/02.Cibersortx/299/CIBERSORTx_299.CD_in.csv", sep=",", head=T)
rownames(CD_in) <- CD_in[,1]
CD_in <- CD_in[,-1]
CD_in <- t(as.matrix(CD_in))
CD_in <- CD_in[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/03.Plot/Cibersortx_299.CD_in.heatmap.png", width=60, height=50, units="cm", res=200)
heatmap <- Heatmap(CD_in, width=unit(40, "cm"), height=unit(40,"cm"), row_names_gp = gpar(fontsize = 30))
draw(heatmap,  heatmap_legend_side="left")
dev.off()


#CD_no
CD_no <- read.table("/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/02.Cibersortx/299/CIBERSORTx_299.CD_no.csv", sep=",", head=T)
rownames(CD_no) <- CD_no[,1]
CD_no <- CD_no[,-1]
CD_no <- t(as.matrix(CD_no))
CD_no <- CD_no[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/03.Plot/Cibersortx_299.CD_no.heatmap.png", width=60, height=50, units="cm", res=200)
heatmap <- Heatmap(CD_no, width=unit(40, "cm"), height=unit(40,"cm"), row_names_gp = gpar(fontsize = 30))
draw(heatmap,  heatmap_legend_side="left")
dev.off()


#UC_in
UC_in <- read.table("/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/02.Cibersortx/299/CIBERSORTx_299.UC_in.csv", sep=",", head=T)
rownames(UC_in) <- UC_in[,1]
UC_in <- UC_in[,-1]
UC_in <- t(as.matrix(UC_in))
UC_in <- UC_in[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/03.Plot/Cibersortx_299.UC_in.heatmap.png", width=60, height=50, units="cm", res=200)
heatmap <- Heatmap(UC_in, width=unit(40, "cm"), height=unit(40,"cm"), row_names_gp = gpar(fontsize = 30))
draw(heatmap,  heatmap_legend_side="left")
dev.off()

#UC_no
UC_no <- read.table("/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/02.Cibersortx/299/CIBERSORTx_299.UC_no.csv", sep=",", head=T)
rownames(UC_no) <- UC_no[,1]
UC_no <- UC_no[,-1]
UC_no <- t(as.matrix(UC_no))
UC_no <- UC_no[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/03.Plot/Cibersortx_299.UC_no.heatmap.png", width=60, height=50, units="cm", res=200)
heatmap <- Heatmap(UC_no, width=unit(40, "cm"), height=unit(40,"cm"), row_names_gp = gpar(fontsize = 30))
draw(heatmap,  heatmap_legend_side="left")
dev.off()

#155
#CD_in
CD_in <- read.table("/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/02.Cibersortx/155/CIBERSORTx_155.CD_in.csv", sep=",", head=T)
rownames(CD_in) <- CD_in[,1]
CD_in <- CD_in[,-1]
CD_in <- t(as.matrix(CD_in))
CD_in <- CD_in[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/03.Plot/Cibersortx_155.CD_in.heatmap.png", width=60, height=50, units="cm", res=200)
heatmap <- Heatmap(CD_in, width=unit(40, "cm"), height=unit(40,"cm"), row_names_gp = gpar(fontsize = 30))
draw(heatmap,  heatmap_legend_side="left")
dev.off()

#CD_no
CD_no <- read.table("/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/02.Cibersortx/155/CIBERSORTx_155.CD_no.csv", sep=",", head=T)
rownames(CD_no) <- CD_no[,1]
CD_no <- CD_no[,-1]
CD_no <- t(as.matrix(CD_no))
CD_no <- CD_no[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/03.Plot/Cibersortx_155.CD_no.heatmap.png", width=60, height=50, units="cm", res=200)
heatmap <- Heatmap(CD_no, width=unit(40, "cm"), height=unit(40,"cm"), row_names_gp = gpar(fontsize = 30))
draw(heatmap,  heatmap_legend_side="left")
dev.off()

#UC_in
UC_in <- read.table("/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/02.Cibersortx/155/CIBERSORTx_155.UC_in.csv", sep=",", head=T)
rownames(UC_in) <- UC_in[,1]
UC_in <- UC_in[,-1]
UC_in <- t(as.matrix(UC_in))
UC_in <- UC_in[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/03.Plot/Cibersortx_155.UC_in.heatmap.png", width=60, height=50, units="cm", res=200)
heatmap <- Heatmap(UC_in, width=unit(40, "cm"), height=unit(40,"cm"), row_names_gp = gpar(fontsize = 30))
draw(heatmap,  heatmap_legend_side="left")
dev.off()

#UC_no
UC_no <- read.table("/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/02.Cibersortx/155/CIBERSORTx_155.UC_no.csv", sep=",", head=T)
rownames(UC_no) <- UC_no[,1]
UC_no <- UC_no[,-1]
UC_no <- t(as.matrix(UC_no))
UC_no <- UC_no[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22),]
png(filename="/scratch/hpc46a05/gy_RNA/IBD/06.Cibersortx/NotfilterByExpr/03.Plot/Cibersortx_155.UC_no.heatmap.png", width=60, height=50, units="cm", res=200)
heatmap <- Heatmap(UC_no, width=unit(40, "cm"), height=unit(40,"cm"), row_names_gp = gpar(fontsize = 30))
draw(heatmap,  heatmap_legend_side="left")
dev.off()
