#!/bin/sh
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/hpc46a05/PBS.OU
#PBS -e /scratch/hpc46a05/Error
#PBS -N hisat2-0
module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate RNAseq

RNAseq="/scratch/hpc46a05/gy_RNA"
REFERENCE="${RNAseq}/Reference"
path="${RNAseq}/IBD"
raw_data_folder="${path}/01.RawData"
hisat2_index_hg38="${REFERENCE}/hisat2_index/hg38"
salmon_index_hg38="${REFERENCE}/salmon_index/hg38"
result_folder="${path}/02.hg38_hisat2_result"
for A in $(cat ${path}/01.Samplelist/split/x0)
do
mkdir ${result_folder}/${A}
hisat2 -q -p 68 -x ${hisat2_index_hg38}/Homo_sapiens.GRCh38.dna.primary_assembly \
    -1 ${raw_data_folder}/${A}_1.fastq.gz \
    -2 ${raw_data_folder}/${A}_2.fastq.gz \
    -S ${result_folder}/${A}/hisat2.dta.sam \
    --downstream-transcriptome-assembly \
    2> ${result_folder}/${A}/hisat2.dta.log

echo $A done
    samtools view -b ${result_folder}/${A}/hisat2.dta.sam \
    | samtools sort -o ${result_folder}/${A}/hisat2.dta.sorted.bam
    samtools index ${result_folder}/${A}/hisat2.dta.sorted.bam
done
