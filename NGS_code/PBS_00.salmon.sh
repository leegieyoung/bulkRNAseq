#!/bin/sh
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/hpc46a05/PBS.OU
#PBS -e /scratch/hpc46a05/Error
#PBS -N salmon-9
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
result_folder="${path}/02.hg38_salmon_result"

mkdir ${path}/02.hg38_salmon_result
for A in $(cat ${path}/01.Samplelist/split/x9)
do
mkdir ${result_folder}/${A}
    salmon quant \
    -i ${salmon_index_hg38} \
    --threads 64 \
    -l A \
    -1 ${raw_data_folder}/${A}_1.fastq.gz \
    -2 ${raw_data_folder}/${A}_2.fastq.gz \
    --validateMappings \
    -o ${result_folder}/${A}/ \
    --gcBias
done
