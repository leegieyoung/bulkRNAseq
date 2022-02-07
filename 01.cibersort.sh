#!/bin/sh
INPUT=$1
mkdir ${INPUT}
CSV_folder="/scratch/hpc46a05/gy_RNA/IBD/04.CSV"
gtf_folder="/scratch/hpc46a05/gy_RNA/Reference/gtf"
summary_folder="${INPUT}"
#step1
sed -n -e '1p' ${INPUT}.csv > raw_${INPUT}.head
echo "gene_id" > head
paste -d ',' head raw_${INPUT}.head > ${INPUT}.head
sed '1d' ${INPUT}.csv > raw_${INPUT}.csv
cat ${INPUT}.head raw_${INPUT}.csv > temp_${INPUT}.csv
sed 's/,/ /g' temp_${INPUT}.csv > ${INPUT}.tab

#step2 annotation of ensembl
awk '{print $1}' ${INPUT}.tab | sed '1d' > ${INPUT}.ensembl
for A in $(cat ${INPUT}.ensembl)
do
grep -w "$A" /scratch/hpc46a05/gy_RNA/Reference/gtf/anno_Homo_sapiens.GRCh38.100.gtf >> ${INPUT}.anno
done
sed -i '1i\ensembl_id gene_symbol function' ${INPUT}.anno

#paste
paste -d ' ' ${INPUT}.anno ${INPUT}.tab > ${INPUT}.result
sed -i 's/ /\t/g' ${INPUT}.result
mv ${INPUT}.result ${summary_folder}/${INPUT}.tab
sed 's/\t\t/\t/g' ${summary_folder}/${INPUT}.tab
sed 's/\t/,/g' ${summary_folder}/${INPUT}.tab > ${summary_folder}/${INPUT}.csv
rm head
rm *.head
rm raw_*
rm temp_*
rm *.anno

#miRNA
sed -n -e '1p' ${summary_folder}/${INPUT}.tab > ${summary_folder}/head
awk '$9 < 0.05' ${summary_folder}/${INPUT}.tab > ${summary_folder}/temp.fdr
awk '$5 > 1 || $5 < -1 {print $0}' ${summary_folder}/temp.fdr_lFC

mkdir ${summary_folder}/01.protein_Coding
grep -w 'protein_coding' ${summary_folder}/temp.fdr_lFC > ${summary_folder}/01.protein_Coding/raw_protein_Coding.tab
cat head ${summary_folder}/01.protein_Coding/raw_protein_Coding.tab > ${summary_folder}/01.protein_Coding/protein_Coding.tab

mkdir ${summary_folder}/02.miRNA
grep -w 'miRNA' ${summary_folder}/temp.fdr_lFC > ${summary_folder}/02.miRNA/raw_miRNA.tab
cat head ${summary_folder}/02.miRNA/raw_miRNA.tab > ${summary_folder}/02.miRNA/miRNA.tab

rm head
rm ${summary_folder}/temp.*
