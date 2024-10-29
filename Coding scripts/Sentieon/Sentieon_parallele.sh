
# DNAseq
# Step 1 Map reads to reference
## Config setup

## sentieon bwa
```
#!/bin/bash

set -e  # 任何命令失败时退出脚本
set -x  # 打印执行的命令

export SENTIEON_LICENSE=/data2/Fulgent_Genetics_cluster_usb361.lic
cd /data2/sentieon-genomics-202112.04/bin

fasta=/data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
readgroup="WGS20220823"
platform="ILLUMINA"
nthreads=10
bam_option="--bam_compression 1"
output_dir=/data1/20240820_NGS_result/WGS_analysis/ 

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
export fasta readgroup platform nthreads bam_option output_dir

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

parallel --link test :::: sample-name :::: fastq1-file :::: fastq2-file :::: output-bam
```

### -j，并行的任务数量
### 这样写的parallel遇到bug，无法正确识别args，!这里应该输入4个数组!
### parallel echo {} ::: arg1 arg2 arg3，每个arg都会被引用到{}的位置并执行
### param1=(p1_1 p1_2 p1_3)
### param2=(p2_1 p2_2 p2_3)
### param3=(p3_1 p3_2 p3_3)
### param4=(p4_1 p4_2 p4_3)
### parallel -j3 fun ::: ${param1[@]} ${param2[@]} ${param3[@]} ${param4[@]}


# Step 2 Calculate data metrics
## Calculate data metrics
#### generate 5 statistical summaries of the data quality and the pipeline data analysis quality results
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
  ./sentieon driver -t 32 -r $fasta -i $bam \
    --algo GCBias --summary "$GC_SUMMARY_TXT" "$GC_METRIC_TXT" \
    --algo MeanQualityByCycle "$MQ_METRIC_TXT" \
    --algo QualDistribution "$QD_METRIC_TXT" \
    --algo InsertSizeMetricAlgo "$IS_METRIC_TXT" \
    --algo AlignmentStat "$ALN_METRIC_TXT"
### --algo refers to algorithms

## generate the plots from the statistical summaries
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

## keep bam and bai files for Sentieon downstream analysis
for bam in "/admin_file/QiaojinLIN/BioCRF/Genomics_Core/202408_NGS_raw/Aligned_data_20240826/WGS_analysis/*.bam";do 
  cp $bam /data1/20240820_NGS_result/WGS_analysis/
done

for bai in "/admin_file/QiaojinLIN/BioCRF/Genomics_Core/202408_NGS_raw/Aligned_data_20240826/WGS_analysis/*.bai";do 
  cp $bai /data1/20240820_NGS_result/WGS_analysis/
done

# Step 3 Remove or mark duplicates
## remove or mark duplicates
```shell
for SORTED_BAM in /data1/20240820_NGS_result/WGS_analysis/*.bam;do
  sample_name=$(basename "$SORTED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$SORTED_BAM")
  ./sentieon driver -t 32 -i "$SORTED_BAM" --algo LocusCollector --fun score_info ${dir_name}/${sample_name}_SCORE.gz
  ./sentieon driver -t 32 -i "$SORTED_BAM" --algo Dedup --rmdup --score_info ${dir_name}/${sample_name}_SCORE.gz --metrics ${dir_name}/${sample_name}_DEDUP_METRIC.txt ${dir_name}/${sample_name}_DEDUPED.bam
done
### output:SCORE.gz;DEDUPED.bam

## Indel realignment (optional)
### 识别可能需要重新比对的区域，这些区域通常是由于存在插入或缺失变异导致的错配，并进行实际的重新比对，以改善对齐质量并减少由于 INDEL 导致的假阳性变异调用
for DEDUPED_BAM in /data1/20240820_NGS_result/WGS_analysis/*DEDUPED.bam;do
  sample_name=$(basename "$DEDUPED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$DEDUPED_BAM")
  ./sentieon driver -t 32 -r $fasta -i $DEDUPED_BAM --algo Realigner ${dir_name}/${sample_name}_realigned.bam
done
### optional:[-k KNOWN_SITES]; the location of the VCF file used as a set of known sites. Multiple collections of known sites can be added by repeating the -k KNOWN_SITES option.
### output: realigned.bam
```
# Step 4 Base quality score recalibration (BQSR)
## Base quality score recalibration (BQSR)
### recalibration math depends on platform (PL) tag of the ReadGroup
### 质量得分（通常以 Phred 得分表示），BQSR 通过分析已知的变异位点来识别和校正bam文件中每个碱基的质量得分中的系统误差
```shell
for REALIGNED_BAM in /data1/20240820_NGS_result/WGS_analysis/*realigned.bam;do
  sample_name=$(basename "$REALIGNED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$REALIGNED_BAM")
  ./sentieon driver -t 32 -r $fasta -i $REALIGNED_BAM --algo QualCal ${dir_name}/${sample_name}_recal_data.txt
done
### optional:[-k KNOWN_SITES]; the location of the VCF file used as a set of known sites. Multiple collections of known sites can be added by repeating the -k KNOWN_SITES option.
### output reacal_data.txt
```

## applies the recalibration to calculate the post calibration data table
## [--algo ReadWriter RECALIBRATED_BAM] additionally apply the recalibration on the BAM file
for REALIGNED_BAM in /data1/20240820_NGS_result/WGS_analysis/*realigned.bam;do
  sample_name=$(basename "$REALIGNED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$REALIGNED_BAM")
  RECAL_TABLE=${dir_name}/${sample_name}_recal_data.txt
  RECAL_TABLE_POST=${dir_name}/${sample_name}_recal_data_post.txt
  ./sentieon driver -t 32 -r $fasta -i $REALIGNED_BAM -q $RECAL_TABLE --algo QualCal $RECAL_TABLE_POST
## creates the data for plotting
  RECAL_RESULT=${dir_name}/${sample_name}_recal_result.csv
  ./sentieon driver -t 32 --algo QualCal --plot --before $RECAL_TABLE --after $RECAL_TABLE_POST $RECAL_RESULT
## plots the calibration data tables
  BQSR_PDF=${dir_name}/${sample_name}_bqsr.pdf
  ./sentieon plot QualCal -o $BQSR_PDF $RECAL_RESULT
done


# Step 5 Variant calling----time consuming
### use GATK Haplotyper to call SNP and indels
for REALIGNED_BAM in /data1/20240820_NGS_result/WGS_analysis/*realigned.bam;do
  sample_name=$(basename "$REALIGNED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$REALIGNED_BAM")
  RECAL_TABLE=${dir_name}/${sample_name}_recal_data.txt
  VARIANT_VCF=${dir_name}/${sample_name}_variant.vcf
  ./sentieon driver -t 32 -r $fasta -i $REALIGNED_BAM -q $RECAL_TABLE  --algo Haplotyper $VARIANT_VCF
done
### dbSNP: the location of the Single Nucleotide Polymorphism database (dbSNP) that will be used to label known variants. Only one dbSNP file.
### What is the difference between genotyper and haplotyper?
### HaplotypeCaller much better at calling indels than position-based callers like UnifiedGenotyper.

# How to perform parallel task? When should we use for loop?
# parallel: 当有多个独立的样本数据需要相同的处理流程时，例如多个样本的质控、比对、排序和去重复等步骤。
# parallel: 当单个样本数据处理可以被分解为多个较小的、可以并行执行的任务时，例如在进行变异检测时，可以将基因组分成多个区域，每个区域由不同的进程或线程处理。
# for:当你需要遍历一个已知大小的数据集合，例如列表、数组或文件中的行。
# for:当你需要对数据集合中的每个元素执行相同操作时，例如对列表中的每个元素进行计算或对文件中的每行进行解析。
# for:当你需要在循环中使用索引时，例如在数组或矩阵计算中。


# 2. DNAscope
# 2.1 Germline variant calling with a machine learning model
## call variants and to apply the machine learning model
### input BAM file should come from a pipeline where only alignment and deduplication have been performed
# PCRFREE=true #PCRFREE=true means the sample is PCRFree, change it to false for PCR samples. In this analysis, 5 PCR loops, non-PCR free library

for DEDUPED_BAM in /data1/20240820_NGS_result/WGS_analysis/*DEDUPED.bam;do
  sample_name=$(basename "$DEDUPED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$DEDUPED_BAM")
  TMP_VARIANT_VCF=${dir_name}/${sample_name}_tmp_variant.vcf
  VARIANT_VCF=${dir_name}/${sample_name}_scope_variant.vcf
  echo "$DEDUPED_BAM $TMP_VARIANT_VCF $VARIANT_VCF" >> /data1/20240820_NGS_result/WGS_analysis/variant_call_tasks.txt
done

# if [ "$PCRFREE" = true ] ; then
# ./sentieon driver -t 32 -r $fasta -i DEDUPED_BAM --algo DNAscope [ -d dbSNP ] --pcr_indel_model none --model ML_MODEL TMP_VARIANT_VCF
# else
# ./sentieon driver -t 32 -r $fasta -i $DEDUPED_BAM --algo DNAscope --model ML_MODEL $TMP_VARIANT_VCF
# fi
# ./sentieon driver -t 32 -r $fasta --algo DNAModelApply --model ML_MODEL -v $TMP_VARIANT_VCF $VARIANT_VCF

vim /data2/sentieon-genomics-202112.04/bin/run_variant_call.sh
#!/bin/bash
tmp_variant_call() {
  local DEDUPED_BAM=$1
  local TMP_VARIANT_VCF=$2
  ./sentieon driver -t 32 -r $fasta -i $DEDUPED_BAM --algo DNAscope --model $ML_MODEL $TMP_VARIANT_VCF
}
variant_call() {
  local TMP_VARIANT_VCF=$2
  local VARIANT_VCF=$3
  ./sentieon driver -t 32 -r $fasta --algo DNAModelApply --model $ML_MODEL -v $TMP_VARIANT_VCF $VARIANT_VCF
}
fasta=/data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
ML_MODEL=/data2/sentieon-genomics-202112.04/bin/sentieon_model/dnascope.model
DEDUPED_BAM="$1"
TMP_VARIANT_VCF="$2"
tmp_variant_call "$DEDUPED_BAM" "$TMP_VARIANT_VCF"
variant_call "$TMP_VARIANT_VCF" "$VARIANT_VCF"

parallel -j 3 ./run_variant_call.sh ::: /data1/20240820_NGS_result/WGS_analysis/variant_call_tasks.txt

### 并行任务还是不行，但这次是函数写错了，运行没有结果，噢，原来是ML_MODEL是个location参数，不是内置的方法，以后每次写复杂任务前都看看核心函数能不能运行
### 问问超哥这个model的文件我们有没有？云上面下载下来的EB的WGS用的ML模型是一个bundle文件，这个要怎么用呢？
### Error: Failed to decode model data
### Finally! It works! No, it doesn't...

### try for loop then
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

## Structural variant calling
### detect >50 nt variant such as insertion, deletion, inversion, translocation
### need to add output break-end (BND) info
### if BQSR is performed, it is possible to input the recalibration table
for DEDUPED_BAM in /data1/20240820_NGS_result/WGS_analysis/*DEDUPED.bam;do
  sample_name=$(basename "$DEDUPED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$DEDUPED_BAM")
  RECAL_TABLE=${dir_name}/${sample_name}_recal_data.txt
  TMP_VARIANT_VCF=${dir_name}/${sample_name}_SV_tmp_variant.vcf
  STRUCTURAL_VARIANT_VCF=${dir_name}/${sample_name}_SV_variant.vcf
  ./sentieon driver -t 32 -r $fasta -i $DEDUPED_BAM -q $RECAL_TABLE --algo DNAscope --var_type bnd $TMP_VARIANT_VCF
  ./sentieon driver -t 32 -r $fasta --algo SVSolver -v $TMP_VARIANT_VCF $STRUCTURAL_VARIANT_VCF
done
### The key difference between a regular VCF and a GVCF is that the GVCF has records for all sites, whether there is a variant call there or not
### The goal is to have every site represented in the file in order to do joint analysis of a cohort in subsequent steps















# 测试一下并行
# With --link you can link the input sources and get one argument from each input source:
# parallel --link echo ::: A B C ::: D E F从变量组中各提取一个，不对等长度变量组合会wrap
# :::传递变量 ::::传递变量文件，:::+会将传递符前后的变量逐个连接，不对等长度变量组合会ignore
# GNU parallel will normally treat a full line as a single argument: It uses \n as argument delimiter. This can be changed with -d:
# parallel -d _ echo :::: abc_-file
# The replacement string {.} removes the extension:
# parallel echo {.} ::: A/B.C
# The replacement string {/} removes the path:
# parallel echo {/} ::: A/B.C
# The replacement string {/.} removes the path and the extension:
# parallel echo {/.} ::: A/B.C
# The replacement string {//} keeps only the path:
# parallel echo {//} ::: A/B.C

#!/bin/bash

set -e  # 任何命令失败时退出脚本
set -x  # 打印执行的命令
for REALIGNED_BAM in /data1/20240820_NGS_result/WGS_analysis/*realigned.bam;do
  sample_name=$(basename "$REALIGNED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$REALIGNED_BAM")
  prefix=${dir_name}/${sample_name}
  echo "$REALIGNED_BAM" >> bam-file
  echo "$prefix" >> prefix-file
done
test() {
  bam=$1
  prefix=$2
  echo "Starting processing for $bam with prefix $prefix"
  ./sentieon driver -t 32 -r /data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -i $bam --algo QualCal ${prefix}_recal_data.txt
  echo "Finished processing for $bam"
}
export -f test
export SENTIEON_LICENSE=/data2/Fulgent_Genetics_cluster_usb361.lic
parallel --link test :::: bam-file :::: prefix-file

### yeah! Finally