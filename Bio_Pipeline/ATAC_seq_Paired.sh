#!/bin/bash
#SBATCH -J Fully_Process_Scripts_ATAC_seq
#SBATCH -o ./output/output-%j.out
#SBATCH -e ./error/error-%j.out
#SBATCH -p standard-s --mem=64G

#Change the directory as necessary
LOG_DIR="/lustre/work/client/users/zeyul/BIT/hg38/data/log/"
HG38_REF_DIR="/lustre/work/client/users/zeyul/tools/hg38_ref_bowtie2/hg38"
DATA_DIR="/lustre/work/client/users/zeyul/BIT/hg38/data/"
FastQC_DIR="/lustre/work/client/users/zeyul/BIT/hg38/FastQC/"
PEAK_DIR="/lustre/work/client/users/zeyul/BIT/hg38/results/"
DiffBind_DIR="/lustre/work/client/users/zeyul/BIT/hg38/DiffBind/"

MIDDLE_DIR="CTCF_KO/"

#Specific the final ATAC-seq filename and SRX number!
ATAC_NAME_1="CTCF_KO_Rep1"
SRX_NUMBER="SRRXXXXXXXX"

mkdir -p "${DATA_DIR}${MIDDLE_DIR}"
mkdir -p "${LOG_DIR}${MIDDLE_DIR}"
mkdir -p "${FastQC_DIR}${MIDDLE_DIR}"
mkdir -p "${PEAK_DIR}${MIDDLE_DIR}"
mkdir -p "${DiffBind_DIR}${MIDDLE_DIR}"

##Download Data
fasterq-dump --progress "$SRX_NUMBER" -o "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.fastq"

##FastQC before Trimming
##FastQC before Trimming
fastqc --noextract --nogroup -o "${FastQC_DIR}${MIDDLE_DIR}" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_1.fastq" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_2.fastq"

##Trimming with trim-galore
trim_galore --paired -q 20 --phred33 --length 25 -e 0.1 --stringency 4 -o "${DATA_DIR}${MIDDLE_DIR}" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_1.fastq" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_2.fastq"

##FastQC post-trimming
fastqc --noextract --nogroup -o "${FastQC_DIR}${MIDDLE_DIR}" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_1_val_1.fq" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_2_val_2.fq"

##bowtie2 Alignment
bowtie2 --very-sensitive -X 2000 -p 8 -q --local \
-x "${HG38_REF_DIR}" \
-1 "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_1_val_1.fq" -2 "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_2_val_2.fq" \
-S "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.sam"

##Samtools transfer to BAM
module load gcc/11.2.0
module load samtools/1.15.1-tyq2r6q

samtools view -h -S -b \
-o "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.bam" \
"${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.sam"

samtools sort -n -o "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_sorted.bam" -O BAM "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.bam"

samtools fixmate -m "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_sorted.bam" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_fixmate.bam"

rm "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_sorted.bam"

samtools sort -o "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_sorted.bam" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_fixmate.bam"

samtools markdup -r -s "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_sorted.bam" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_Final.bam"

samtools index "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_Final.bam"

rm "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_1.fastq"
rm "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_2.fastq"
rm "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.sam"
rm "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_sorted.bam"
rm "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_fixmate.bam"

macs2 callpeak -t "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_Final.bam" -f BAM -g mm -n "${PEAK_DIR}${MIDDLE_DIR}${ATAC_NAME_1}" -B -q 0.01 


#For DiffBind
mv "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_Final.bam" "${DiffBind_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_Final.bam"
mv "${PEAK_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_peaks.narrowPeak" "${DiffBind_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_peaks.narrowPeak"




