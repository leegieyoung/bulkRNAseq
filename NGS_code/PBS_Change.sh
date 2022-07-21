#!/bin/sh
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/hpc46a05/PBS/PBS.OU/
#PBS -e /scratch/hpc46a05/PBS/Error/
#PBS -N Change
module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate RNAseq

Input=$1
RNAseq="/scratch/hpc46a05/gy_RNA"
REFERENCE="${RNAseq}/Reference"
raw_data_folder="/scratch/hpc46a05/raw_data/01.RNAseq/02.Colon_Cancer_raw.data/fastq"
hisat2_index_hg38="${REFERENCE}/hisat2_index/hg38"
salmon_index_hg38="${REFERENCE}/salmon_index/hg38"
result_folder="${RNAseq}/colon_Cancer/01.hisat/"

mkdir ${result_folder}/Change
hisat2 -q -p 68 -x ${hisat2_index_hg38}/Homo_sapiens.GRCh38.dna.primary_assembly \
    -1 ${raw_data_folder}/Change.R1.fq.gz \
    -2 ${raw_data_folder}/Change.R2.fq.gz \
    -S ${result_folder}/Change/hisat2.dta.sam \
    --downstream-transcriptome-assembly \
    2> ${result_folder}/Change/hisat2.dta.log

samtools view -b ${result_folder}/Change/hisat2.dta.sam | samtools sort -o ${result_folder}/Change/hisat2.dta.sorted.bam
samtools index ${result_folder}/Change/hisat2.dta.sorted.bam

module purge
module load singularity/3.9.7
singularity exec /scratch/hpc46a05/image/RNAseq.sif Rscript /scratch/hpc46a05/gy_RNA/Code/hg38/hisat2/PBS_Colon_Cancer/R/PBS_Change.R 
