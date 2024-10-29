### **Date**
>20240924~20240925
### **Environment**
>server:10.30.62.1
### **Task**
* AVITI RNA-seq fastq alignment
* RNA-seq downstream analysis pipeline test
* 10X pipeline test
* ONT pipeline test


### **Log**
#### AVITI RNA-seq fastq alignment

为测序样本选择正确的参考基因组（GRCm39）进行alignment

```shell
cd /data1/Analysis/AVITI_analysis_20240923
# mouse sample(YF3)
nextflow run nf-core/rnaseq \
     -profile docker \
     --input samplesheet.csv \
     --outdir output \
     --gtf /data1/Data/Reference/Mus_musculus.GRCm39.112.gtf.gz \
     --fasta /data1/Data/Reference/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz \
     -c nf_custom.config
# human sample(Li)
nextflow run nf-core/rnaseq \
     -profile docker \
     --input samplesheet_Li_human.csv \
     --outdir output/Li_human \
     --gtf /data1/Data/Reference/Homo_sapiens.GRCh38.112.gtf.gz \
     --fasta /data1/Data/Reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz \
     -c nf_custom.config
```

#### Pipeline test
##### Dataset download
```shell
# 下载10X的reference
cd /data1/Data/Reference/Reference_10X
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
# 10X multiplex example dataset(use cellranger multi pipeline)
cd /data1/Data/ExampleData
curl -O "wget https://cg.10xgenomics.com/samples/cell-exp/6.0.0/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Multiplex/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Multiplex_fastqs.tar"
tar -xf SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Multiplex_fastqs.tar
# 10X singleplex example dataset(use cellranger count pipeline)
curl -O "https://cf.10xgenomics.com/samples/cell-exp/6.1.0/3p_Citrate_CPT/3p_Citrate_CPT_fastqs.tar"
tar -xf 3p_Citrate_CPT_fastqs.tar
```
##### Run pipeline
```cellranger multi``` pipeline inputs:
* ```--id``` name the output directory
* ```--csv``` points to the FASTQ files, and contains other parameters
```shell
cd /data1/Data/ExampleData/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K
mkdir inputs
cd inputs
nano SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K.csv
# copy paste
[gene-expression]
ref,/data1/Data/Reference/Reference_10X/refdata-gex-GRCh38-2024-A
create-bam,true

[libraries]
fastq_id,fastqs,lanes,feature_types
SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_1_gex,/data1/Data/ExampleData/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_1_gex,any,Gene Expression
SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_1_multiplexing_capture,/data1/Data/ExampleData/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_1_multiplexing_capture,any,Multiplexing Capture

[samples]
sample_id,cmo_ids,description
Jurkat,CMO301,Jurkat
Raji,CMO302,Raji
# ctrl+O to save file and ctrl+X exit nano editor.

# Check cellranger multi parameters use:
cellranger multi-template --parameters
cellranger multi --help

# Run cellranger multi
cellranger multi --id=Jurkat_Raji_10K --csv=SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K.csv
# starts at 17:41:21, finished at 18:37:18, output save at /data1/Data/ExampleData/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/inputs/Jurkat_Raji_10K/outs/
# Count matrix output saved at /data1/Data/ExampleData/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/inputs/Jurkat_Raji_10K/outs/multi/count
# Summary saved at /data1/Data/ExampleData/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/inputs/Jurkat_Raji_10K/outs/per_sample_outs/
```


(测试)SRA数据(GEO)fetch
```shell
# (test) Fetch SRA files
cd /data1/Data/Reference/
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
export PATH=$PATH:$PWD/sratoolkit.3.1.1-ubuntu64/bin
prefetch srafile
```
* Note
  
[10X产品相关名词解释](https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-glossary)


