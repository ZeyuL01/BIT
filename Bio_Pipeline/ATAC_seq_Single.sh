#!/bin/bash
#SBATCH -J Fully_Process_Scripts_ATAC_seq
#SBATCH -o ./output/output-%j.out
#SBATCH -e ./error/error-%j.out
#SBATCH -p standard-s --mem=64G

#Change the directory as necessary
LOG_DIR="/lustre/work/client/users/zeyul/BIT/mm10/data/log/"
MM10_REF_DIR="/lustre/work/client/users/zeyul/tools/mm10_ref_bowtie2/mm10"
DATA_DIR="/lustre/work/client/users/zeyul/BIT/mm10/data/"
FastQC_DIR="/lustre/work/client/users/zeyul/BIT/mm10/FastQC/"
PEAK_DIR="/lustre/work/client/users/zeyul/BIT/mm10/results/"
DiffBind_DIR="/lustre/work/client/users/zeyul/BIT/mm10/DiffBind/"

MIDDLE_DIR="RUNX1_KO/"

#Specific the final ATAC-seq filename and SRX number!
ATAC_NAME_1="RUNX1_KO_Rep1"
SRX_NUMBER="SRRXXXXXX"

mkdir -p "${DATA_DIR}${MIDDLE_DIR}"
mkdir -p "${LOG_DIR}${MIDDLE_DIR}"
mkdir -p "${FastQC_DIR}${MIDDLE_DIR}"
mkdir -p "${PEAK_DIR}${MIDDLE_DIR}"
mkdir -p "${DiffBind_DIR}${MIDDLE_DIR}"

##Download Data
fasterq-dump --progress "$SRX_NUMBER" -o "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.fastq"

##FastQC before Trimming
fastqc --noextract --nogroup -o "${FastQC_DIR}${MIDDLE_DIR}" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.fastq"

##Trimming with trim-galore
trim_galore -q 20 --phred33 --length 25 -e 0.1 --stringency 4 -o "${DATA_DIR}${MIDDLE_DIR}" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.fastq"

##FastQC post-trimming
fastqc --noextract --nogroup -o "${FastQC_DIR}${MIDDLE_DIR}" "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_trimmed.fq"

##bowtie2 Alignment
bowtie2 --very-sensitive -X 2000 -p 8 -q --local \
-x "${MM10_REF_DIR}" \
-U "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_trimmed.fq" \
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

rm "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.fastq"
rm "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}.sam"
rm "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_sorted.bam"
rm "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_fixmate.bam"

macs2 callpeak -t "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_Final.bam" -f BAM -g mm -n "${PEAK_DIR}${MIDDLE_DIR}${ATAC_NAME_1}" -B -q 0.01

#For DiffBind
mv "${DATA_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_Final.bam" "${DiffBind_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_Final.bam"
mv "${PEAK_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_peaks.narrowPeak" "${DiffBind_DIR}${MIDDLE_DIR}${ATAC_NAME_1}_peaks.narrowPeak"
