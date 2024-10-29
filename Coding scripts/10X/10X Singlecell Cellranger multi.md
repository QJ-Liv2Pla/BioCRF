### （Example）Dataset download
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
## Step 1 Fastq to count matrix
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
# starts at 17:41:21, finished at 18:37:18

# To explore output
cd Jurkat_Raji_10K/outs
tree

# Count matrix output saved at /data1/Data/ExampleData/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/inputs/Jurkat_Raji_10K/outs/multi/count
# Summary saved at /data1/Data/ExampleData/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/inputs/Jurkat_Raji_10K/outs/per_sample_outs/. Warnings will show if results are abnormal.

```

## Step 3 Single-cell downstream analysis
#### In R
```R
library(Seurat)
library(tidyverse)
library(ggplot2)

# 1. Read 10X data
## Read in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix.
My10Xproject <- Read10X(data.dir = "/data1/Data/ExampleData/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/inputs/Jurkat_Raji_10K/outs/per_sample_outs/Raji/count/sample_filtered_feature_bc_matrix")
Sobj <- CreateSeuratObject(counts = My10Xproject, project = "Jurkat_Raji_10Xmulti", min.cells = 3, min.features = 200)

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
```

```R
## Remove low quality cells, 
Sobj <- subset(Sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

```R
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
```

```R
# 5. Cells clustering
## Use Louvain algorithm to iteratively group cells together
Sobj <- FindNeighbors(Sobj, dims = 1:10)
## Decide the granularity of cell clustering. Higher resolution value leads to more coarse clustering result.
Sobj <- FindClusters(Sobj, resolution = 0.5)
## Visualization of clustering result based on UMAP algorithm
Sobj <- RunUMAP(Sobj, dims = 1:10)
DimPlot(Sobj, reduction = "umap")
```
```R
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
```

```R
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

#### Last edited by Qijian on 20240927