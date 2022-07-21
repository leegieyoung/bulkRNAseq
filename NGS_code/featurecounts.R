library(Rsubread)

gtf="/scratch/hpc46a05/gy_RNA/Reference/gtf/Homo_sapiens.GRCh38.100.gtf.gz"
pre_dir="/scratch/hpc46a05/gy_RNA/colon_Cancer/01.hisat/"


FeatureCounts <-function(bam){
	bamfile <-paste0(pre_dir,bam,"/","hisat2.dta.sorted.bam") 
	fct=featureCounts(file=bamfile, isPairedEnd=TRUE, countMultiMappingReads=TRUE, countReadPairs=TRUE, isGTFAnnotationFile=TRUE, nthreads=68, annot.ext=gtf, tmpDir="/scratch/hpc46a05/gy_RNA/IBD/02.featureCounts_result/temp")
	write.table(fct$counts, file=paste0("/scratch/hpc46a05/gy_RNA/colon_Cancer/02.featureCounts_result/",bam,".txt"), row.names=T, col.names=T, sep=' ', quote=F)
}
