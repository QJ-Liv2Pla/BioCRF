# 10X Singlecell Analysis
## Step 1 Download and install Cellranger and Seurat

```shell
# Download and install Cellranger
cd /opt
sudo mkdir Cellranger
curl -o cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1725283779&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=AM7BT31MtszJgCXdi2VtjcpMTlsdoEV5sLSNFj9he7uztN9k~A-5CTQbNGa9U8i4Y373t8HqtVguCKHwfe8tifFPnn3ZtJXm24VDXIsgriotwA4ypTvD5xUZGtfhax1bQv1wFzhBtoCqsrHkQFjpoJB6FgxVOn1n6UO4zT3E3LRov-OanjF5qoGooaCHXBGE95xGsl4lT7nub8rapAzP5~K0Jp3FcNoiCYsx6YMdUeJy18989zRc~X7fOZYenfwgIBwAcDe0rRkYypQR59asxQqHZudBgOHFuUWZGbHlVyvD7JlCSs99ji9yz86bnRGN7qHcTkaD6gw8NNC0paWizQ__"
tar -xzvf cellranger-x.y.z.tar.gz
export PATH=/opt/Cellranger/cellranger-8.0.1:$PATH

# Download human reference genome
cd /data1/Data/Reference
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"

# Create conda env for Singlecell analysis and install R and Seurat
source /opt/anaconda/anaconda3/bin/activate
conda create -n Singlecell
conda activate Singlecell
conda install conda-forge::r-base --yes
R
```
```R
# In R
install.packages('Seurat')
```

## Step 2 Cellranger analysis
### (Optional) If raw files are SRA files in .sra format, extracts data in FASTQ-format
```shell
SRAdir=$"/Path/to/srafile/dir"
Outdir=$"/output/dir/"
fasterq-dump --split-files $SRAdir --include-technical -O $Outdir
# This step generates 3 fastq files: (1)Sample index read，(2)Read1 and (3)read2. 

# Rename them as followed:
Sample_ID_I1_001.fastq
Sample_ID_R1_001.fastq
Sample_ID_R2_001.fastq

# For cellranger analysis, zip them as followed:
for file in /home/qijian/scRNAseq/A8_PRJNA634159/fastq/SRR11821876/*.fastq;do pigz $file;done
```
### Use cellranger to generate QC and count matrix from FASTQ
```shell
# Get example data
cd /data1/Data/ExampleData
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar

# Run cellranger count
cellranger count --id=run_count_1kpbmcs \
   --fastqs=/data1/Data/ExampleData/pbmc_1k_v3_fastqs \
   --sample=pbmc_1k_v3 \
   --transcriptome=/data1/Data/Reference/refdata-gex-GRCh38-2024-A

# For multiple fastqs
Project_id=$"Projectid"
Reference=$"/Path/to/GRCh38/ref/file"
Outdir=$"/Path/to/outputdir/for/countmatrix"
for i in /Directory-to-fastq/*; do
     cellranger count --id=$Project_id \
          --transcriptome=$Reference \
          --fastqs=$i \
          --create-bam=true \
          --localcores=64 \
          --localmem=512 \
          --output-dir=$Outdir;
done
```

## Step 3 Single-cell downstream analysis
```R
library(Seurat)
library(tidyverse)
library(ggplot2)

# 1. Read 10X data
## Read in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix.
My10Xproject <- Read10X(data.dir = "/Path/To/10X/countmatrix/")
Sobj <- CreateSeuratObject(counts = My10Xproject, project = "My10Xproject", min.cells = 3, min.features = 200)

## Read in metadata which contains cell index and other info.
metadata <- "/Path/to/metadata.csv"
MyMetaData <- read.csv(metadata, sep = ",", header = TRUE, fill = TRUE)
Sobj@meta.data <- MyMetaData

# 2. Pre-processing
## QC visualization
Sobj[["percent.mt"]] <- PercentageFeatureSet(Sobj, pattern = "^MT-")
VlnPlot(Sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(Sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
## Remove low quality cells, 
Sobj <- subset(Sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 3. Data normalization
## Normalize based on the assumption that each cell originally contains the same number of RNA molecules. Otherwise use SCTransform()
Sobj <- NormalizeData(Sobj, normalization.method = "LogNormalize", scale.factor = 10000)
Sobj <- NormalizeData(Sobj)
## Feature selection to reduce data dimensionality
Sobj <- FindVariableFeatures(Sobj, selection.method = "vst", nfeatures = 2000)
Sobj <- ScaleData(Sobj, features = VariableFeatures(Sobj))

# 4. Dimension reduction
Sobj <- RunPCA(Sobj, features = VariableFeatures(object = Sobj), npcs = 50)
DimPlot(Sobj, reduction = "pca")
## Chose number of PCs for downstream analysis based on major standadr deviation change in elbow plot
ElbowPlot(Sobj)

# 5. Cells clustering
## Use Louvain algorithm to iteratively group cells together
Sobj <- FindNeighbors(Sobj, dims = 1:10)
## Decide the granularity of cell clustering. Higher resolution value leads to more coarse clustering result.
Sobj <- FindClusters(Sobj, resolution = 0.5)
## Visualization of clustering result based on UMAP algorithm
Sobj <- RunUMAP(Sobj, dims = 1:10)
DimPlot(Sobj, reduction = "umap")

# 6. DEG analysis
## Analyze the differential expressed genes in all clusters
Sobj.markers <- FindAllMarkers(Sobj, only.pos = TRUE)
Sobj.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
VlnPlot(Sobj, features = top10)
FeaturePlot(Sobj, features = Sobj, features = top10)

# 7. Assign cell names to cluster
new.cluster.ids <- c("T cells", "B cells", "NK cells", "DC")
names(new.cluster.ids) <- levels(Sobj)
Sobj <- RenameIdents(Sobj, new.cluster.ids)
## Plot
plot <- DimPlot(Sobj, reduction = "umap", label = TRUE, label.size = 4.5) + 
     xlab("UMAP 1") + 
     ylab("UMAP 2") +
     theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + 
     guides(colour = guide_legend(override.aes = list(size = 10)))
## Save result
ggsave(filename = "../output/images/MyProject_CellClustering.jpg", height = 7, width = 12, plot = plot, quality = 50)
```

#### Last edited by Qijian on 20240920




