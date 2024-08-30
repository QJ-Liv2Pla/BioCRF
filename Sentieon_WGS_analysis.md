# DNAseq
# Step 1 Map reads to reference
## Config setup
```shell
# omit
```
## sentieon bwa

```shell
#!/bin/bash
export SENTIEON_LICENSE=/data2/Fulgent_Genetics_cluster_usb361.lic
cd /data2/sentieon-genomics-202112.04/bin

set -e  # exit script when encouter failure
set -x  # print command

fasta=/data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
readgroup="WGS20220823"
platform="ILLUMINA"
nthreads=10
bam_option="--bam_compression 1"
output_dir=/data1/20240820_NGS_result/WGS_analysis/ 

export fasta readgroup platform nthreads bam_option output_dir

run_bam() {
  sample_name=$1
  fastq_1=$2
  fastq_2=$3
  output_bam=$4
  echo "Starting processing for $sample_name, read fastq from $fastq_1 and $fastq_2 "
  (./sentieon bwa mem -M -R "@RG\tID:$readgroup\tSM:$sample_name\tPL:$platform" -t $nthreads -K 10000000 $fasta $fastq_1 $fastq_2 || echo -n 'error') | ./sentieon util sort $bam_option -r $fasta -o $output_bam -t $nthreads --sam2bam -i -
  echo "Finished processing for $sample_name, output bam file $output_bam"
}

export -f run_bam

for dir in /admin_file/BioCRF/Genomics_Core/Element_AVITI/202408_NGS_raw/WGS/WGS_merge/*/; do
  sample_name=$(basename $dir)
  fastq_1="${dir}/${sample_name}_R1.fastq.gz"
  fastq_2="${dir}/${sample_name}_R2.fastq.gz"
  output_bam="${output_dir}${sample_name}_sorted.bam"
  echo "$sample_name" >> sample-name
  echo "$fastq_1" >> fastq1-file
  echo "$fastq_2" >> fastq2-file
  echo "$output_bam" >> output-bam
done

parallel --link run_bam :::: sample-name :::: fastq1-file :::: fastq2-file :::: output-bam
```
#### output: .bam & .bam.bai

# Step 2 Calculate data metrics
## Calculate data metrics
```shell
for bam in /data1/20240820_NGS_result/WGS_analysis/*.bam; do
  sample_name=$(basename "$bam" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$bam")
  GC_SUMMARY_TXT="${dir_name}/${sample_name}_gc_summary.txt"
  GC_METRIC_TXT="${dir_name}/${sample_name}_gc_metric.txt"
  MQ_METRIC_TXT="${dir_name}/${sample_name}_mq_metric.txt"
  QD_METRIC_TXT="${dir_name}/${sample_name}_qd_metric.txt"
  IS_METRIC_TXT="${dir_name}/${sample_name}_is_metric.txt"
  ALN_METRIC_TXT="${dir_name}/${sample_name}_aln_metric.txt"
# generate 5 statistical summaries of the data quality and the pipeline data analysis quality results
  ./sentieon driver -t 32 -r $fasta -i $bam \
    --algo GCBias --summary "$GC_SUMMARY_TXT" "$GC_METRIC_TXT" \
    --algo MeanQualityByCycle "$MQ_METRIC_TXT" \
    --algo QualDistribution "$QD_METRIC_TXT" \
    --algo InsertSizeMetricAlgo "$IS_METRIC_TXT" \
    --algo AlignmentStat "$ALN_METRIC_TXT"
# generate the plots from the statistical summaries
  GC_METRIC_PDF="${dir_name}/${sample_name}_gc_metric.pdf"
  MQ_METRIC_PDF="${dir_name}/${sample_name}_mq_metric.pdf"
  QD_METRIC_PDF="${dir_name}/${sample_name}_qd_metric.pdf"
  IS_METRIC_PDF="${dir_name}/${sample_name}_is_metric.pdf"
  ./sentieon plot GCBias -o "$GC_METRIC_PDF" "$GC_METRIC_TXT"
  ./sentieon plot MeanQualityByCycle -o "$MQ_METRIC_PDF" "$MQ_METRIC_TXT"
  ./sentieon plot QualDistribution -o "$QD_METRIC_PDF" "$QD_METRIC_TXT"
  ./sentieon plot InsertSizeMetricAlgo -o "$IS_METRIC_PDF" "$IS_METRIC_TXT"
done
```
#### output: GC_SUMMARY_TXT, GC_METRIC_TXT, MQ_METRIC_TXT, QD_METRIC_TXT, IS_METRIC_TXT, ALN_METRIC_TXT, etc.
# Step 3 Remove or mark duplicates
## remove or mark duplicates
```shell
for SORTED_BAM in /data1/20240820_NGS_result/WGS_analysis/*.bam;do
  sample_name=$(basename "$SORTED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$SORTED_BAM")
  ./sentieon driver -t 32 -i "$SORTED_BAM" --algo LocusCollector --fun score_info ${dir_name}/${sample_name}_SCORE.gz
  ./sentieon driver -t 32 -i "$SORTED_BAM" --algo Dedup --rmdup --score_info ${dir_name}/${sample_name}_SCORE.gz --metrics ${dir_name}/${sample_name}_DEDUP_METRIC.txt ${dir_name}/${sample_name}_DEDUPED.bam
done
```
#### output: SCORE.gz;DEDUPED.bam

## Indel realignment (optional)
### 识别可能需要重新比对的区域，这些区域通常是由于存在插入或缺失变异导致的错配，并进行实际的重新比对，以改善对齐质量并减少由于INDEL导致的假阳性变异调用
```shell
for DEDUPED_BAM in /data1/20240820_NGS_result/WGS_analysis/*DEDUPED.bam;do
  sample_name=$(basename "$DEDUPED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$DEDUPED_BAM")
  ./sentieon driver -t 32 -r $fasta -i $DEDUPED_BAM --algo Realigner ${dir_name}/${sample_name}_realigned.bam
done
# optional:[-k KNOWN_SITES]; the location of the VCF file used as a set of known sites. Multiple collections of known sites can be added by repeating the -k KNOWN_SITES option.
```
#### output: realigned.bam
# Step 4 Base quality score recalibration (BQSR)
### 质量得分（通常以Phred得分表示），BQSR 通过分析已知的变异位点来识别和校正bam文件中每个碱基的质量得分中的系统误差
```shell
for REALIGNED_BAM in /data1/20240820_NGS_result/WGS_analysis/*realigned.bam;do
  sample_name=$(basename "$REALIGNED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$REALIGNED_BAM")
  ./sentieon driver -t 32 -r $fasta -i $REALIGNED_BAM --algo QualCal ${dir_name}/${sample_name}_recal_data.txt
done
```
#### output: reacal_data.txt

# Step 5 Variant calling
### use GATK Haplotyper to call SNP and indels
```shell

for REALIGNED_BAM in /data1/20240820_NGS_result/WGS_analysis/*realigned.bam;do
  sample_name=$(basename "$REALIGNED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$REALIGNED_BAM")
  RECAL_TABLE=${dir_name}/${sample_name}_recal_data.txt
  VARIANT_VCF=${dir_name}/${sample_name}_variant.vcf
  ./sentieon driver -t 32 -r $fasta -i $REALIGNED_BAM -q $RECAL_TABLE  --algo Haplotyper $VARIANT_VCF
done
```
#### (1) --dbSNP: the location of the Single Nucleotide Polymorphism database (dbSNP) that will be used to label known variants. Only one dbSNP file.
#### (2) HaplotypeCaller is much better at calling indels than position-based callers like UnifiedGenotyper.
#### output: recall.txt, variant.vcf

# DNAScope
# Step 1 (After DNAseq)Germline variant calling with a machine learning model
### Input BAM file should come from a pipeline where only alignment and deduplication have been performed. If PCRFree samples are used, add "--pcr_indel_model".
```shell
fasta=/data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
ML_MODEL=/data2/sentieon-genomics-202112.04/bin/sentieon_model/dnascope.model
for DEDUPED_BAM in /data1/20240820_NGS_result/WGS_analysis/*DEDUPED.bam;do
  sample_name=$(basename "$DEDUPED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$DEDUPED_BAM")
  TMP_VARIANT_VCF=${dir_name}/${sample_name}_tmp_variant.vcf
  VARIANT_VCF=${dir_name}/${sample_name}_scope_variant.vcf
  ./sentieon driver -t 32 -r $fasta -i $DEDUPED_BAM --algo DNAscope --model $ML_MODEL $TMP_VARIANT_VCF
  ./sentieon driver -t 32 -r $fasta --algo DNAModelApply --model $ML_MODEL -v $TMP_VARIANT_VCF $VARIANT_VCF
done
```
#### output: scope_variant.vcf

# Step 2 Structural variant calling
### SV is used to detect >50 nt variant such as insertion, deletion, inversion, translocation. Break-end (BND) info is needed. If BQSR is performed, it is possible to input the recalibration table
```shell
for DEDUPED_BAM in /data1/20240820_NGS_result/WGS_analysis/*DEDUPED.bam;do
  sample_name=$(basename "$DEDUPED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$DEDUPED_BAM")
  RECAL_TABLE=${dir_name}/${sample_name}_recal_data.txt
  TMP_VARIANT_VCF=${dir_name}/${sample_name}_SV_tmp_variant.vcf
  STRUCTURAL_VARIANT_VCF=${dir_name}/${sample_name}_SV_variant.vcf
  ./sentieon driver -t 32 -r $fasta -i $DEDUPED_BAM -q $RECAL_TABLE --algo DNAscope --var_type bnd $TMP_VARIANT_VCF
  ./sentieon driver -t 32 -r $fasta --algo SVSolver -v $TMP_VARIANT_VCF $STRUCTURAL_VARIANT_VCF
done
```
#### The key difference between a regular VCF and a GVCF is that the GVCF has records for all sites, whether there is a variant call there or not. The goal is to have every site represented in the file in order to do joint analysis of a cohort in subsequent steps
#### output: recal_data.txt, SV_variant.vcf

