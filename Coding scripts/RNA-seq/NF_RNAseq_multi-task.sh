# Nextflow RNA-seq analysis
# Step 1
### Create a samplesheet in csv format which contain information about sample name, file dir of fastq files and strandness of library.
vim /data1/20240820_NGS_result/RNA_seq_analysis/20240826_input/samplesheet.csv
sample,fastq_1,fastq_2,strandedness
RNA_32_1,/admin_file/BioCRF/Genomics_Core/Element_AVITI/202408_NGS_raw/RNA/RNA_32_1/RNA_32_1_R1.fastq.gz,/admin_file/BioCRF/Genomics_Core/Element_AVITI/202408_NGS_raw/RNA/RNA_32_1/RNA_32_1_R2.fastq.gz,auto
RNA_32_2,/admin_file/BioCRF/Genomics_Core/Element_AVITI/202408_NGS_raw/RNA/RNA_32_2/RNA_32_2_R1.fastq.gz,/admin_file/BioCRF/Genomics_Core/Element_AVITI/202408_NGS_raw/RNA/RNA_32_2/RNA_32_2_R2.fastq.gz,auto

# Step 2
### Customize the configuration of analysis. For detailed syntax, visit:
[configuration syntax](https://www.nextflow.io/docs/latest/config.html#process-selectors)
vim /data1/20240820_NGS_result/RNA_seq_analysis/20240826_input/custom_20240826.config
executor {
     name = 'local'
     maxParallel = 20
}
process {
     cpus = 64
     memory = 500.GB
}
# Step 3
nextflow run nf-core/rnaseq \
     -profile docker \
     --input /data1/20240820_NGS_result/RNA_seq_analysis/20240826_input/samplesheet.csv \
     --outdir /data1/20240820_NGS_result/RNA_seq_analysis/20240826_output \
     --gtf /data1/Data/Reference/Homo_sapiens.GRCh38.112.gtf.gz \
     --fasta /data1/Data/Reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
     -c /data1/20240820_NGS_result/RNA_seq_analysis/20240826_input/custom_20240826.config
## Note:
### (1) Be careful wheather adapter were trimmed before NF workflows
### (2) Do not use -c <file> to specify parameters as this will result in errors. 
### (3) Custom config files specified with -c must only be used for tuning process resource specifications, other infrastructural tweaks (such as output directories), or module arguments (args).

# Pipeline settings can be provided in a yaml or json file via -params-file <file>.
# with params.yaml containing:
input: <SAMPLESHEET>
outdir: <OUTDIR>
genome: 'GRCh37'
<...>

### pipeline will create the following files in working directory:
#### work  Directory containing the nextflow working files
#### <OUTDIR>  Finished results in specified location (defined with --outdir)
#### .nextflow_log  Log file from Nextflow
#### Other nextflow hidden files, eg. history of pipeline runs and old logs.
#### Last edited by qijian 20240920
