#!/bin/bash

# ===============================================================
# Germline Variant Calling Pipeline (GATK4 Best Practices)
# Whole Genome Sequencing (Paired-End 2×100 bp)
# Author: Adithya
# NOTE: Demonstration pipeline 
# ===============================================================

set -e

# -----------------------------
# Directory Setup
# -----------------------------
mkdir -p reads support_files/hg38 ngs_variant_project/{reads,aligned_reads,results,data}

# -----------------------------
# STEP 0: Download FASTQ Reads
# -----------------------------
echo "Downloading FASTQ files..."
wget -P reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

# ===============================================================
# PREP FILES (Run Once)
# ===============================================================

# -----------------------------
# Download Reference Genome (hg38)
# -----------------------------
echo "Downloading hg38 reference..."
wget -P support_files/hg38 https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip support_files/hg38/hg38.fa.gz

# Index reference (.fai)
samtools faidx support_files/hg38/hg38.fa

# Dictionary (.dict)
gatk CreateSequenceDictionary -R support_files/hg38/hg38.fa -O support_files/hg38/hg38.dict

# Known sites (dbSNP for BQSR)
echo "Downloading known-sites for BQSR..."
wget -P support_files/hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P support_files/hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# ===============================================================
# STEP 1: FASTQC
# ===============================================================
echo "STEP 1: QC with FastQC"
fastqc reads/SRR062634_1.filt.fastq.gz -o ngs_variant_project/reads/
fastqc reads/SRR062634_2.filt.fastq.gz -o ngs_variant_project/reads/

# ===============================================================
# STEP 2: BWA-MEM ALIGNMENT
# ===============================================================
echo "STEP 2: Mapping with BWA-MEM"

bwa index support_files/hg38/hg38.fa

bwa mem -t 4 \
  -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" \
  support_files/hg38/hg38.fa \
  reads/SRR062634_1.filt.fastq.gz \
  reads/SRR062634_2.filt.fastq.gz \
  > ngs_variant_project/aligned_reads/SRR062634.paired.sam

# ===============================================================
# STEP 3: MARK DUPLICATES & SORT
# ===============================================================
echo "STEP 3: MarkDuplicates & Sorting"

gatk MarkDuplicatesSpark \
  -I ngs_variant_project/aligned_reads/SRR062634.paired.sam \
  -O ngs_variant_project/aligned_reads/SRR062634_sorted_dedup_reads.bam

# ===============================================================
# STEP 4: BASE QUALITY SCORE RECALIBRATION (BQSR)
# ===============================================================
echo "STEP 4: Base Quality Recalibration"

gatk BaseRecalibrator \
  -I ngs_variant_project/aligned_reads/SRR062634_sorted_dedup_reads.bam \
  -R support_files/hg38/hg38.fa \
  --known-sites support_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
  -O ngs_variant_project/data/recal_data.table

gatk ApplyBQSR \
  -I ngs_variant_project/aligned_reads/SRR062634_sorted_dedup_reads.bam \
  -R support_files/hg38/hg38.fa \
  --bqsr-recal-file ngs_variant_project/data/recal_data.table \
  -O ngs_variant_project/aligned_reads/SRR062634_sorted_dedup_bqsr_reads.bam

# ===============================================================
# STEP 5: ALIGNMENT METRICS
# ===============================================================
echo "STEP 5: Alignment & Insert Metrics"

gatk CollectAlignmentSummaryMetrics \
  R=support_files/hg38/hg38.fa \
  I=ngs_variant_project/aligned_reads/SRR062634_sorted_dedup_bqsr_reads.bam \
  O=ngs_variant_project/aligned_reads/alignment_metrics.txt

gatk CollectInsertSizeMetrics \
  INPUT=ngs_variant_project/aligned_reads/SRR062634_sorted_dedup_bqsr_reads.bam \
  OUTPUT=ngs_variant_project/aligned_reads/insert_size_metrics.txt \
  HISTOGRAM_FILE=ngs_variant_project/aligned_reads/insert_size_histogram.pdf

# ===============================================================
# STEP 6: VARIANT CALLING (HaplotypeCaller)
# ===============================================================
echo "STEP 6: Variant Calling"

gatk HaplotypeCaller \
  -R support_files/hg38/hg38.fa \
  -I ngs_variant_project/aligned_reads/SRR062634_sorted_dedup_bqsr_reads.bam \
  -O ngs_variant_project/results/raw_variants.vcf

# Extract SNPs and Indels
gatk SelectVariants -R support_files/hg38/hg38.fa -V ngs_variant_project/results/raw_variants.vcf --select-type SNP -O ngs_variant_project/results/raw_snps.vcf
gatk SelectVariants -R support_files/hg38/hg38.fa -V ngs_variant_project/results/raw_variants.vcf --select-type INDEL -O ngs_variant_project/results/raw_indels.vcf


echo "Variant Calling Sucessfull"