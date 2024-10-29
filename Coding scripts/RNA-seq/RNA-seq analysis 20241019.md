# RNA-seq analysis 20241019


# 1. Output
## 1.1 Bases2Fastq Project QC Report
/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/20240927_Wu_Group_QC.html
## 1.2
### 1.2.1 MultiQC
/data1/Analysis/NGS_analysis_20241019/output/multiqc/star_salmon/multiqc_report.html
### 1.2.2 RNA quantification
/data1/Analysis/NGS_analysis_20241019/output/star_salmon


# 2. Manuscript
mkdir /data1/Analysis/NGS_analysis_20241019
cd /data1/Analysis/NGS_analysis_20241019
```shell
sample,fastq_1,fastq_2,strandedness
C_1_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_1_1/C_1_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_1_1/C_1_1_R2.fastq.gz,auto
C_2_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_2_1/C_2_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_2_1/C_2_1_R2.fastq.gz,auto
C_3_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_3_1/C_3_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_3_1/C_3_1_R2.fastq.gz,auto
L_1_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_1_1/L_1_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_1_1/L_1_1_R2.fastq.gz,auto
L_2_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_2_1/L_2_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_2_1/L_2_1_R2.fastq.gz,auto
L_3_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_3_1/L_3_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_3_1/L_3_1_R2.fastq.gz,auto
N_1_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_1_1/N_1_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_1_1/N_1_1_R2.fastq.gz,auto
N_2_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_2_1/N_2_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_2_1/N_2_1_R2.fastq.gz,auto
N_3_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_3_1/N_3_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_3_1/N_3_1_R2.fastq.gz,auto
N_4_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_4_1/N_4_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_4_1/N_4_1_R2.fastq.gz,auto

vim /data1/Analysis/NGS_analysis_20241019/custom.config
executor {   
     name = 'local'
     maxParallel = 20
}
process {
     cpus = 64
     memory = 500.GB
}
```

```shell
nextflow run nf-core/rnaseq \
     -profile docker \
     --input samplesheet.csv \
     --outdir output \
     --gtf /data1/Data/Reference/Mus_musculus.GRCm39.112.gtf.gz \
     --fasta /data1/Data/Reference/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz \
     -c custom.config
```

  
```shell
cd /data1/Analysis/NGS_analysis_20241019/abundance_analysis
vim samplesheet.csv
sample,fastq_1,fastq_2,condition,replicate
C_1_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_1_1/C_1_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_1_1/C_1_1_R2.fastq.gz,C,1
C_2_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_2_1/C_2_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_2_1/C_2_1_R2.fastq.gz,C,2
C_3_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_3_1/C_3_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/C_3_1/C_3_1_R2.fastq.gz,C,3
L_1_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_1_1/L_1_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_1_1/L_1_1_R2.fastq.gz,L,1
L_2_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_2_1/L_2_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_2_1/L_2_1_R2.fastq.gz,L,2
L_3_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_3_1/L_3_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/L_3_1/L_3_1_R2.fastq.gz,L,3
N_1_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_1_1/N_1_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_1_1/N_1_1_R2.fastq.gz,N,1
N_2_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_2_1/N_2_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_2_1/N_2_1_R2.fastq.gz,N,2
N_3_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_3_1/N_3_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_3_1/N_3_1_R2.fastq.gz,N,3
N_4_1,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_4_1/N_4_1_R1.fastq.gz,/staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/N_4_1/N_4_1_R2.fastq.gz,N,4

vim contrast.csv
id,variable,reference,target
C_vs_L,condition,C,L
C_vs_N,condition,C,N
N_vs_L,condition,N,L

nextflow run nf-core/differentialabundance \
     -profile docker \
     --input /data1/Analysis/NGS_analysis_20241019/abundance_analysis/samplesheet.csv \
     --contrasts /data1/Analysis/NGS_analysis_20241019/abundance_analysis/contrasts.csv \
     --matrix /data1/Analysis/NGS_analysis_20241019/output/star_salmon/salmon.merged.gene_counts.tsv \
     --transcript_length_matrix /data1/Analysis/NGS_analysis_20241019/output/star_salmon/salmon.merged.gene_lengths.tsv \
     --gtf /data1/Data/Reference/Mus_musculus.GRCm39.112.gtf.gz \
     --outdir /data1/Analysis/NGS_analysis_20241019/abundance_analysis/output3 \
     --gprofiler2_run true \
     --gprofiler2_organism mmusculus \
     --gene_sets_files /data1/Analysis/NGS_analysis_20241019/abundance_analysis/gprofiler_full_mmusculus.name.gmt
```


```shell
sudo find /staging/data/fastq/20241016_AV242212_20241016_B/Samples/20240927_Wu_Group/* -type d -exec cp -r {}   \
```