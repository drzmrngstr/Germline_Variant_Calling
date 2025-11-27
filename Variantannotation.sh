# ===============================================================
# STEP 7: VARIANT FILTERING (Hard-filtering as demonstration)
# ===============================================================
# NOTE: For WGS it's recommended to use VQSR when you have many samples.
# Here we apply GATK hard filters suitable for a single-sample demo.

echo "STEP 7: Variant Filtering (SNPs and INDELs)"

# Filter SNPs
gawk -v OFS='	' 'BEGIN{print "# SNP FILTERING"}' > /dev/null

gatk VariantFiltration \
  -R support_files/hg38/hg38.fa \
  -V ngs_variant_project/results/raw_snps.vcf \
  -O ngs_variant_project/results/filtered_snps.vcf \
  --filter-name "QD_filter" --filter-expression "QD < 2.0" \
  --filter-name "FS_filter" --filter-expression "FS > 60.0" \
  --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
  --filter-name "SOR_filter" --filter-expression "SOR > 4.0" \
  --filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
  --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
  --genotype-filter-expression "DP < 10" --genotype-filter-name "DP_filter" \
  --genotype-filter-expression "GQ < 10" --genotype-filter-name "GQ_filter"

# Filter INDELs

gatk VariantFiltration \
  -R support_files/hg38/hg38.fa \
  -V ngs_variant_project/results/raw_indels.vcf \
  -O ngs_variant_project/results/filtered_indels.vcf \
  --filter-name "QD_filter" --filter-expression "QD < 2.0" \
  --filter-name "FS_filter" --filter-expression "FS > 200.0" \
  --filter-name "SOR_filter" --filter-expression "SOR > 10.0" \
  --genotype-filter-expression "DP < 10" --genotype-filter-name "DP_filter" \
  --genotype-filter-expression "GQ < 10" --genotype-filter-name "GQ_filter"

# Select variants that PASS filters

echo "Selecting PASS variants"

gatk SelectVariants \
  --exclude-filtered \
  -V ngs_variant_project/results/filtered_snps.vcf \
  -O ngs_variant_project/results/passed-snps.vcf

gatk SelectVariants \
  --exclude-filtered \
  -V ngs_variant_project/results/filtered_indels.vcf \
  -O ngs_variant_project/results/passed-indels.vcf

# excluding those variants that failed the genotype filters

echo "Removing variants that failed genotype filters"


cat ngs_variant_project/results/passed-snps.vcf|grep -v -E "DP_filter|GQ_filter" >ngs_variant_project/results/passed-snps-gtfiltered.vcf
cat ngs_variant_project/results/passed-indels.vcf| grep -v -E "DP_filter|GQ_filter" > ngs_variant_project/results/passed-indels-gtfiltered.vcf

# ===============================================================
# STEP 8: ANNOTATION (Funcotator)
# ==================S=============================================
echo "STEP 8: Annotation with Funcotator (optional)"

# Download funcotator data sources (one-time; large). Uncomment if needed.
#download and copy it to the gatk folder
#By default funcotator comes with ClinVar and gencode annotation if you want gnoMAD untar the gnoMAD file present.


./gatk FuncotatorDataSourceDownloader --germline --hg38 --validate-integrity --extract-after-download

gatk Funcotator \
--variant ngs_variant_project/results/passed-snps-gtfiltered.vcf \
--reference support_files/hg38/hg38.fa \
--ref-version hg38 \
--data-sources-path /home/mrrngstr/apps/gatk-4.6.2.0/funcotator_dataSources.v1.8.hg38.20230908g \
--output ngs_variant_project/results/passed-snps-gtfiltered-functotated.vcf \
--output-file-format VCF

gatk Funcotator \
--variant ngs_variant_project/results/passed-indels-gtfiltered.vcf \
--reference support_files/hg38/hg38.fa  \
--ref-version hg38 \
--data-sources-path /home/mrrngstr/apps/gatk-4.6.2.0/funcotator_dataSources.v1.8.hg38.20230908g \
--output ngs_variant_project/results/passed-indels-gtfiltered-functotated.vcf \
--output-file-format VCF

# ===============================================================
# STEP 9: Convert vcf to  table
# ===============================================================
echo "STEP 9: Export variants to table"

# If Funcotator output exists, extract FUNCOTATION field; otherwise use basic fields from VCF
if [ -f ngs_variant_project/results/passed-snps-gtfiltered-functotated.vcf ]; then
  gatk VariantsToTable \
    -V ngs_variant_project/results/passed-snps-gtfiltered-functotated.vcf \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O ngs_variant_project/results/output_snps.table
else
  gatk VariantsToTable \
    -V ngs_variant_project/results/passed-snps-gtfiltered.vcf \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AN -F DP -F AF \
    -O ngs_variant_project/results/output_snps.table
fi

if [ -f ngs_variant_project/results/passed-indels-gtfiltered-functotated.vcf ]; then
  gatk VariantsToTable \
    -V ngs_variant_project/results/passed-indels-gtfiltered-functotated.vcf \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O ngs_variant_project/results/output_indels.table
else
  gatk VariantsToTable \
    -V ngs_variant_project/results/passed-indels-gtfiltered.vcf \
    -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F AC -F AN -F DP -F AF \
    -O ngs_variant_project/results/output_indels.table
fi
#The table need further refining which can be done using excel or the terminal itself varints can be sorted and analysed based on gene id and other parameters
# End of Script
echo "Pipeline Completed Successfully."