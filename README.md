GERMLINE VARIANT CALLING USING GATK GOODPRACTICE

📌 Overview

This repository contains a complete Germline Variant Calling Pipeline following GATK4 Best Practices, designed for Whole Genome Sequencing (WGS) paired-end Illumina reads (2×100 bp).
It includes two major components:

variantcalling.sh – Performs QC, alignment, duplicate marking, BQSR, and variant calling.

Variantannotation.sh – Performs SNP/INDEL filtering, annotation with Funcotator, and export to tables.

The pipeline is suitable for demonstration / academic use, running on Linux with commonly used bioinformatics tools (BWA, Samtools, GATK4, FastQC).

🧬 Pipeline Structure
├── variantcalling.sh         # Steps 0–6 (QC → Alignment → BQSR → Variant calling)
├── Variantannotation.sh      # Steps 7–9 (Filtering → Annotation → Export tables)
├── ngs_variant_project/
│   ├── reads/
│   ├── aligned_reads/
│   ├── results/
│   └── data/
└── support_files/
    └── hg38/                 # Reference genome + known sites

⚙️ Requirements

Ensure the following tools are installed and available in your PATH:

GATK4

BWA-MEM

Samtools

FastQC

gawk

wget

Java 8+

🚀 Usage
1. Make scripts executable
chmod +x variantcalling.sh
chmod +x Variantannotation.sh

2. Run variant calling

This script performs:

Downloading FASTQ files

Downloading hg38 reference

FastQC

BWA MEM alignment

Duplicate marking

Base recalibration

Variant calling (HaplotypeCaller)

Splitting raw variants into SNPs and INDELs

./variantcalling.sh

3. Run filtering + annotation

This script performs:

Hard filtering for SNPs/INDELs

Filtering genotypes (DP/GQ filters)

Annotation using Funcotator (optional)

Exporting final VCFs to tabular files

./Variantannotation.sh

📝 Steps Included
🧪 Step 1: Quality Control

FastQC on raw FASTQ files.

🔧 Step 2: Alignment

BWA-MEM indexing + alignment with read groups.

🧹 Step 3: Mark Duplicates

Using MarkDuplicatesSpark.

📉 Step 4: BQSR

BaseRecalibrator and ApplyBQSR using dbSNP known sites.

📊 Step 5: Metrics

Alignment metric summary

Insert size metrics (+ histogram)

🧬 Step 6: Variant Calling

Run HaplotypeCaller

Split into SNP and INDEL sets

🔬 Step 7: Hard Filtering

Includes filters on:

QD

FS

MQ

SOR

MQRankSum

ReadPosRankSum

Genotype DP

Genotype GQ

🧾 Step 8: Annotation (Optional)

Funcotator using GATK dataSources

Generates annotated VCF

📑 Step 9: Export to Tables

Generates tabular output for SNPs and INDELs

Includes FUNCOTATION field when available

📂 Output Structure
Variantcalling:
ngs_variant_project/
  ├── reads/                          # FASTQC output
  ├── aligned_reads/
  │   ├── sorted BAM
  │   ├── metrics.txt
  │   └── insert_size_histogram.pdf
  └── results/
      ├── raw_variants.vcf
      ├── raw_snps.vcf
      └── raw_indels.vcf

Variantannotation:
ngs_variant_project/results/
  ├── filtered_snps.vcf
  ├── filtered_indels.vcf
  ├── passed-snps-gtfiltered.vcf
  ├── passed-indels-gtfiltered.vcf
  ├── *functotated.vcf (optional)
  ├── output_snps.table
  └── output_indels.table

⚠️ Notes & Limitations

This is a demonstration pipeline.
For production WGS pipelines, GATK recommends VQSR instead of hard filters (requires ≥30 samples).

Ensure sufficient disk space (reference genome + data).

Funcotator data sources must be downloaded manually (large files).

🧑‍💻 Author

Adithya (drzmrngstr)
Bioinformatics | Computational Biology | NGS Data Analysis

⭐ If you use this pipeline, consider giving the repository a star!

# Germline_Variant_Calling
