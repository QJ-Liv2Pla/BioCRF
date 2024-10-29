# Nextflow RNA-seq analysis
# Step 1
### Create a samplesheet in csv format which contain information about sample name, file dir of fastq files and strandness of library.
```shell
vim /data1/20240820_NGS_result/RNA_seq_analysis/20240826_input/samplesheet.csv
sample,fastq_1,fastq_2,strandedness
RNA_32_1,/admin_file/BioCRF/Genomics_Core/Element_AVITI/202408_NGS_raw/RNA/RNA_32_1/RNA_32_1_R1.fastq.gz,/admin_file/BioCRF/Genomics_Core/Element_AVITI/202408_NGS_raw/RNA/RNA_32_1/RNA_32_1_R2.fastq.gz,auto
RNA_32_2,/admin_file/BioCRF/Genomics_Core/Element_AVITI/202408_NGS_raw/RNA/RNA_32_2/RNA_32_2_R1.fastq.gz,/admin_file/BioCRF/Genomics_Core/Element_AVITI/202408_NGS_raw/RNA/RNA_32_2/RNA_32_2_R2.fastq.gz,auto
```

# Step 2
### Customize the configuration of analysis. For detailed syntax, visit:
[configuration syntax](https://www.nextflow.io/docs/latest/config.html#process-selectors)
```shell
vim /data1/20240820_NGS_result/RNA_seq_analysis/20240826_input/custom_20240826.config
executor {
     name = 'local'
     maxParallel = 20
}
process {
     cpus = 64
     memory = 500.GB
}
```
# Step 3
```shell
nextflow run nf-core/rnaseq \
     -profile docker \
     --input /data1/20240820_NGS_result/RNA_seq_analysis/20240826_input/samplesheet.csv \
     --outdir /data1/20240820_NGS_result/RNA_seq_analysis/20240826_output \
     --gtf /data1/Data/Reference/Homo_sapiens.GRCh38.112.gtf.gz \
     --fasta /data1/Data/Reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
     -c /data1/20240820_NGS_result/RNA_seq_analysis/20240826_input/custom_20240826.config
```
## Note:
#### (1) Be careful wheather adapter were trimmed before NF workflows
#### (2) Do not use -c <file> to specify parameters as this will result in errors. 
#### (3) Custom config files specified with -c must only be used for tuning process resource specifications, other infrastructural tweaks (such as output directories), or module arguments (args).

## Pipeline settings 
#### (1) Can be provided in a yaml or json file via -params-file <file>. With params.yaml containing:
```shell
input: <SAMPLESHEET>
outdir: <OUTDIR>
genome: 'GRCh37'
<...>
```

#### (2) Pipeline will create the following files in working directory:
```shell
<work directory> # Directory containing the nextflow working files
<OUTDIR>  # Finished results in specified location (defined with --outdir)
.nextflow_log  # Log file from Nextflow
Other nextflow hidden files # eg. history of pipeline runs and old logs.
```

## differentialabundance
#### samplesheet.e.g.
sample,fastq_1,fastq_2,condition,replicate,batch
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,control,1,A
CONTROL_REP2,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,control,2,B
CONTROL_REP3,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,control,3,A
TREATED_REP1,AEG588A2_S1_L002_R1_001.fastq.gz,AEG588A2_S1_L002_R2_001.fastq.gz,treated,1,B
TREATED_REP2,AEG588A2_S1_L003_R1_001.fastq.gz,AEG588A2_S1_L003_R2_001.fastq.gz,treated,2,A
TREATED_REP3,AEG588A2_S1_L004_R1_001.fastq.gz,AEG588A2_S1_L004_R2_001.fastq.gz,treated,3,B
#### contrasts
id,variable,reference,target
Condition_Temp40_Temp37,condition,Temp40,Temp37
Condition_Temp40_Temp32,condition,Temp40,Temp32
Condition_Temp37_Temp32,condition,Temp37,Temp32

```shell
nextflow run nf-core/differentialabundance \
     -profile docker \
     --input samplesheet.csv \
     --contrasts contrasts.csv \
     --matrix /admin_file/QiaojinLIN/BioCRF/Genomics_Core/202408_NGS_raw/Aligned_data_20240826/RNAseq_analysis/star_salmon/salmon.merged.gene_counts.tsv \
     --transcript_length_matrix /admin_file/QiaojinLIN/BioCRF/Genomics_Core/202408_NGS_raw/Aligned_data_20240826/RNAseq_analysis/star_salmon/salmon.merged.gene_lengths.tsv \
     --gtf /data1/Data/Reference/Homo_sapiens.GRCh38.112.gtf.gz \
     --gsea_run true \
     --gene_sets_files /data1/Analysis/NGS_analysis_20240820/DifferentialAbundance/gprofiler_full_hsapiens.name.gmt \
     --outdir /data1/Analysis/NGS_analysis_20240820/DifferentialAbundance/output2
```

```shell
nextflow run nf-core/differentialabundance \
     -profile docker \
     --input samplesheet.csv \
     --contrasts contrasts.csv \
     --matrix /admin_file/QiaojinLIN/BioCRF/Genomics_Core/202408_NGS_raw/Aligned_data_20240826/RNAseq_analysis/star_salmon/salmon.merged.gene_counts.tsv \
     --transcript_length_matrix /admin_file/QiaojinLIN/BioCRF/Genomics_Core/202408_NGS_raw/Aligned_data_20240826/RNAseq_analysis/star_salmon/salmon.merged.gene_lengths.tsv \
     --gtf /data1/Data/Reference/Homo_sapiens.GRCh38.112.gtf.gz \
     --outdir /data1/Analysis/NGS_analysis_20240820/DifferentialAbundance/output3 \
     --gprofiler2_run true \
     --gprofiler2_organism hsapiens \
     --gene_sets_files /data1/Analysis/NGS_analysis_20240820/DifferentialAbundance/gprofiler_full_hsapiens.name.gmt # Necessary
```

##### Last edited by qijian 20241015
