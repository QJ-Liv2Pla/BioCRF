> Be aware of extracting sample name and file dir correctly
> sample info: LinGroup, HelaCells, TwoSamples, ~60GB each
```shell
mkdir /data1/Analysis/Sentieon_DNA_analysis_20241023
cd /data1/Analysis/Sentieon_DNA_analysis_20241023
vim sentieon_20241023_script.sh
```

```shell
#!/bin/bash
# export SENTIEON_LICENSE=ISoGeneticServer-013:8990
export SENTIEON_LICENSE=/data2/Fulgent_Genetics_cluster_usb361.lic
set -e  # exit script when encouter failure
set -x  # print command
fasta=/data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
readgroup="S32C20241023"
platform="ELEMENT"
nthreads=16
bam_option="--bam_compression 1"
output_dir=/data1/Analysis/Sentieon_DNA_analysis_20241023/ # Don't miss the last slash /

export fasta readgroup platform nthreads bam_option output_dir

run_bam() {
  sample_name=$1
  fastq_1=$2
  fastq_2=$3
  output_bam=$4
  echo "Starting processing for $sample_name, read fastq from $fastq_1 and $fastq_2 "
  (/data2/sentieon-genomics-202112.04/bin/./sentieon bwa mem -M -R "@RG\tID:$readgroup\tSM:$sample_name\tPL:$platform" -t $nthreads -K 10000000 $fasta $fastq_1 $fastq_2 || echo -n 'error') | /data2/sentieon-genomics-202112.04/bin/./sentieon util sort $bam_option -r $fasta -o $output_bam -t $nthreads --sam2bam -i -
  echo "Finished processing for $sample_name, output bam file $output_bam"
}

export -f run_bam

for dir in /staging/data/fastq/20241016_AV242212_20241016_B/Samples/202410_Lin_Group/*/; do
  sample_name=$(basename $dir)
  fastq_1="${dir}${sample_name}_R1.fastq.gz"
  fastq_2="${dir}${sample_name}_R2.fastq.gz"
  output_bam="${output_dir}${sample_name}_sorted.bam" # 
  echo "$sample_name" >> sample-name
  echo "$fastq_1" >> fastq1-file
  echo "$fastq_2" >> fastq2-file
  echo "$output_bam" >> output-bam
done

parallel --link run_bam :::: sample-name :::: fastq1-file :::: fastq2-file :::: output-bam
```


```shell
# Error
cmdline: /data2/sentieon-genomics-202112.04/bin/../libexec/util sort --bam_compression 1 -r /data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -o /data1/Analysis/Sentieon_DNA_analysis_2024102332C_3_sorted.bam -t 16 --sam2bam -i -
Failed to contact the license server at ISoGeneticServer-013:8990
Failed to contact the license server at ISoGeneticServer-013:8990
Starting processing for 32C_2, read fastq from /staging/data/fastq/20241016_AV242212_20241016_B/Samples/202410_Lin_Group/32C_2/32C_2_R1.fastq.gz and /staging/data/fastq/20241016_AV242212_20241016_B/Samples/202410_Lin_Group/32C_2/32C_2_R2.fastq.gz
Finished processing for 32C_2, output bam file /data1/Analysis/Sentieon_DNA_analysis_2024102332C_2_sorted.bam
cmdline: /data2/sentieon-genomics-202112.04/bin/../libexec/util sort --bam_compression 1 -r /data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -o /data1/Analysis/Sentieon_DNA_analysis_2024102332C_2_sorted.bam -t 16 --sam2bam -i -
Failed to contact the license server at ISoGeneticServer-013:8990
Failed to contact the license server at ISoGeneticServer-013:8990
```
```shell
./sentieon licsrvr --start --log licsrvr.log -R /data2/sentieon-genomics-202112.04/report.txt /data2/Fulgent_Genetics_cluster_usb361.lic
```

```shell
fasta=/data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
for bam in /data1/Analysis/Sentieon_DNA_analysis_20241023/*.bam; do
  sample_name=$(basename "$bam" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$bam")
  GC_SUMMARY_TXT="${dir_name}/${sample_name}_gc_summary.txt"
  GC_METRIC_TXT="${dir_name}/${sample_name}_gc_metric.txt"
  MQ_METRIC_TXT="${dir_name}/${sample_name}_mq_metric.txt"
  QD_METRIC_TXT="${dir_name}/${sample_name}_qd_metric.txt"
  IS_METRIC_TXT="${dir_name}/${sample_name}_is_metric.txt"
  ALN_METRIC_TXT="${dir_name}/${sample_name}_aln_metric.txt"
# generate 5 statistical summaries of the data quality and the pipeline data analysis quality results
  sentieon driver -t 32 -r $fasta -i $bam \
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
  sentieon plot GCBias -o "$GC_METRIC_PDF" "$GC_METRIC_TXT"
  sentieon plot MeanQualityByCycle -o "$MQ_METRIC_PDF" "$MQ_METRIC_TXT"
  sentieon plot QualDistribution -o "$QD_METRIC_PDF" "$QD_METRIC_TXT"
  sentieon plot InsertSizeMetricAlgo -o "$IS_METRIC_PDF" "$IS_METRIC_TXT"
# remove or mark duplicates
  sentieon driver -t 32 -i $bam --algo LocusCollector --fun score_info ${dir_name}/${sample_name}_SCORE.gz
  sentieon driver -t 32 -i $bam --algo Dedup --rmdup --score_info ${dir_name}/${sample_name}_SCORE.gz --metrics ${dir_name}/${sample_name}_DEDUP_METRIC.txt ${dir_name}/${sample_name}_DEDUPED.bam
done
```

```shell
# Indel realignment (optional)
for DEDUPED_BAM in /data1/Analysis/Sentieon_DNA_analysis_20241023/*DEDUPED.bam;do
  sample_name=$(basename "$DEDUPED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$DEDUPED_BAM")
  sentieon driver -t 32 -r $fasta -i $DEDUPED_BAM --algo Realigner ${dir_name}/${sample_name}_realigned.bam
done

# Base quality score recalibration
for REALIGNED_BAM in /data1/Analysis/Sentieon_DNA_analysis_20241023/*realigned.bam;do
  sample_name=$(basename "$REALIGNED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$REALIGNED_BAM")
  sentieon driver -t 32 -r $fasta -i $REALIGNED_BAM --algo QualCal ${dir_name}/${sample_name}_recal_data.txt
### QualCal: Missing VCF file of known polymorphic sites, this may lead to inaccurate results
# use GATK Haplotyper to call SNP and indels
  RECAL_TABLE=${dir_name}/${sample_name}_recal_data.txt
  VARIANT_VCF=${dir_name}/${sample_name}_variant.vcf
  sentieon driver -t 32 -r $fasta -i $REALIGNED_BAM -q $RECAL_TABLE  --algo Haplotyper $VARIANT_VCF
done

# Germline variant calling with a machine learning model
ML_MODEL=/data2/sentieon-genomics-202112.04/bin/sentieon_model/dnascope.model
for DEDUPED_BAM in /data1/Analysis/Sentieon_DNA_analysis_20241023/*DEDUPED.bam;do
  sample_name=$(basename "$DEDUPED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$DEDUPED_BAM")
  TMP_VARIANT_VCF=${dir_name}/${sample_name}_tmp_variant.vcf
  VARIANT_VCF=${dir_name}/${sample_name}_scope_variant.vcf
  sentieon driver -t 32 -r $fasta -i $DEDUPED_BAM --algo DNAscope --model $ML_MODEL $TMP_VARIANT_VCF
  sentieon driver -t 32 -r $fasta --algo DNAModelApply --model $ML_MODEL -v $TMP_VARIANT_VCF $VARIANT_VCF
done

# Structural variant calling
for DEDUPED_BAM in /data1/Analysis/Sentieon_DNA_analysis_20241023/*DEDUPED.bam;do
  sample_name=$(basename "$DEDUPED_BAM" .bam | awk -F'_' '{print $1"_"$2}')
  dir_name=$(dirname "$DEDUPED_BAM")
  RECAL_TABLE=${dir_name}/${sample_name}_recal_data.txt
  TMP_VARIANT_VCF=${dir_name}/${sample_name}_SV_tmp_variant.vcf
  STRUCTURAL_VARIANT_VCF=${dir_name}/${sample_name}_SV_variant.vcf
  sentieon driver -t 32 -r $fasta -i $DEDUPED_BAM -q $RECAL_TABLE --algo DNAscope --var_type bnd $TMP_VARIANT_VCF
  sentieon driver -t 32 -r $fasta --algo SVSolver -v $TMP_VARIANT_VCF $STRUCTURAL_VARIANT_VCF
done
```


## Use DISCVRSeq/VariantQC to check vcf file.
### 20241029
```shell
docker pull ghcr.io/bimberlab/discvrseq:latest
# Check if docker works
docker run ghcr.io/bimberlab/discvrseq VariantQC --help
fasta=/data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
export fasta
java -jar DISCVRSeq.jar VariantQC \
  -R $fasta \
  -V SimpleExample.vcf.gz \
  -O output.html
```

