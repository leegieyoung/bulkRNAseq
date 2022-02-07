#!/bin/sh
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/hpc46a05/PBS.OU
#PBS -e /scratch/hpc46a05/Error
#PBS -N KB-026-1-TR.sorted.bam

module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate RNAseq

gtf="/scratch/hpc46a05/gy_RNA/Reference/gtf/Homo_sapiens.GRCh38.100.gtf.gz"
temp_folder="/scratch/hpc46a05/gy_RNA/IBD/02.featureCounts_result/temp/"
OUT_folder="/scratch/hpc46a05/gy_RNA/IBD/02.featureCounts_result/"
path="/scratch/hpc46a05/raw_data/01.RNAseq/01.IBD_raw.data/IBD.bam/"

featureCounts ${path}/KB-026-1-TR.sorted.bam -Q 20 -T 64 -a ${gtf} -p -g 'gene_id' --tmpDir ${temp_folder} -o ${OUT_folder}/KB-026-1-TR.sorted.bam.fc_P
