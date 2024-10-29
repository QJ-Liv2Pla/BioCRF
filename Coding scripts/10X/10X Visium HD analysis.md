# Space Ranger

#### Space Ranger is a suite of analysis tools designed for 10x Genomics Visium data. It enables users to profile the whole transcriptome in formalin-fixed & paraffin-embedded (FFPE), fixed frozen (FxF), and fresh frozen (FF) tissues. Space Ranger is compatible with both brightfield and fluorescence microscopy images. Pipelines include:
#### (1) The **spaceranger mkfastq** demultiplexes the Illumina sequencer's base call files (BCLs) for each flow cell directory into FASTQ files;
#### (2) The **spaceranger count** pipeline inputs a reference, a microscope slide image, and FASTQ files to generate feature-barcode matrices, identify clusters, and perform differential gene expression;
#### (3) The **spaceranger aggr** pipeline can be used to aggregate samples into a single feature-barcode matrix.

# Step 1 Download and install Space Ranger
```shell
cd /opt
sudo mkdir Spaceranger
cd Spaceranger
curl -o spaceranger-3.1.1.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-3.1.1.tar.gz?Expires=1726856991&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=PMsrdGplIptYKbsXiHmE-6I6rhVxYe2pH4o-2QbWW05guFO-T3SQJ0B7NTlFEdjPATigONyZvaJhuLTUPIZ1R6JVWY05NB8fAEd2DGS4reHsHIKobmR3XsAvIuUBSZtgw91v0242rCFuFKDS7AIwpujBxHPWHEbF11m23d1TWeKuPIcyqs71rA0PT~mRgnK2NyjSdYszW-p-2ajFjZx-1DFMpHzzL-t35uOn-HQHXtBNL35XBN4GJe-sheIzNig5kiQQDRXkJMUk9bf0mud7u8ZG4mAx1RFyjtIhN~gUf16PuSpt2RNm3TXOSwk7kpc5S5Rtaco68R0wn~kyPghCuA__"
tar -zxvf spaceranger-3.1.0.tar.gz
export PATH=/opt/Spaceranger/spaceranger-3.1.0:$PATH
```

## (Optional Step) Bcl To FASTQ
### Skip this step if raw files are FASTQ

```shell
# spaceranger mkfastq
# 1. Create a CSV sample sheet
vim spaceranger-bcl-sample.csv
Lane,Sample,Index
1,test_sample,SI-TT-D9

# 2. Run mkfastq
Bclpath=$"/PATH/TO/tiny-bcl"
Samplesheet=/PATH/TO/spaceranger-bcl-sample.csv
spaceranger mkfastq --id=tiny-bcl \
                    --run=$Bclpath \
                    --csv=$Samplesheet
```

# Step 2 Visium HD Analysis with spaceranger count
### Following inputs are required:
* **CytAssist image** in ```TIFF``` format (```--cytaimage```)
* **Microscope image** (optional) in either ```TIFF```, ```QPTIFF```, ```BTF```, or ```JPEG``` format:
    * ```--image``` for a brightfield microscope image
    * ```--darkimage``` for a dark background fluorescence microscope image
    * ```--colorizedimage``` for a composite colored fluorescence microscope image
* **Slide parameters** specified by ```--slide``` & ```--area```
* **Reference transcriptome** (```--transcriptome```)
* **probe set CSV** (```--probe-set```)   
#### For ***multiomic*** experiments, separate libraries for the GEX and PEX reads are generated. A ```CSV``` file indicating the input data folder, sample name, and library type of each input library, then pass this file to ```spaceranger count``` using the ```--libraries``` flag.
```shell
Reference=/PATH/TO/REF
Fastqs=/data1/Data/ExampleData/Visium_HD_10X/Visium_HD_Human_Breast_Cancer_Fresh_Frozen_fastqs
ProbeSet=/path/to/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv
CytaImage=/path/to/TIF
Image=/path/to/BTF
Outdir=/data1/Analysis/Visium_HD_analysis_20240929
spaceranger count --id=visium_hd_count \
     --transcriptome=$Reference \
     --fastqs=$Fastqs \
     --probe-set=$ProbeSet \
     --slide=H1-YD7CDZK \
     --area=A1 \
     --cytaimage=$CytaImage \
     --image=$Image \
     --create-bam=false
     --output-dir=Outdir
```
## Step 3 Downstream analysis
### Check spaceranger output
[Spaceranger Output Details](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview)

### Analysis, visualization, and integration of Visium HD spatial dataset
#### 1. Load Visisum HD data
```R
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
```

```R
# The Seurat can store multiple binnings/resolutions in different assays
# bin.size parameter specifies resolutions to load (8 and 16um are loaded by default)
# Switch between resolutions by changing the assay

localdir <- ""
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))

# Setting default assay changes between 8um and 16um binning
Assays(object)
DefaultAssay(object) <- "Spatial.008um"
```

```R
vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot
```

#### 2. Normalization
```R
# Default log-normalization, optimal methods should be chosen
DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object)
```

#### 3. Visualize gene expression
```R
# switch spatial resolution to 16um from 8um
DefaultAssay(object) <- "Spatial.016um"
p1 <- SpatialFeaturePlot(object, features = "Rorb") + ggtitle("Rorb expression (16um)")
```

#### 4. Unsupervised clustering
```R
# Sketch the Visium HD dataset, perform clustering on the subsampled cells, and then project the cluster labels back to the full dataset
# note that data is already normalized
DefaultAssay(object) <- "Spatial.008um"
object <- FindVariableFeatures(object)
object <- ScaleData(object)
# we select 50,0000 cells and create a new 'sketch' assay
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
```

```R
# switch analysis to sketched cells
DefaultAssay(object) <- "sketch"

# perform clustering workflow
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 3)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
```

```R
object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)
```

```R
DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2
```
```R
SpatialDimPlot(object, label = T, repel = T, label.size = 4)
```
```R
Idents(object) <- "seurat_cluster.projected"
cells <- CellsByIdentities(object, idents = c(0, 4, 32, 34, 35))
p <- SpatialDimPlot(object,
  cells.highlight = cells[setdiff(names(cells), "NA")],
  cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p
```
```R
# Crete downsampled object to make visualization either
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p
```

#### 5. Identifying spatially-defined tissue domains
```R
if (!requireNamespace("Banksy", quietly = TRUE)) {
  remotes::install_github("prabhakarlab/Banksy@devel")
}
library(SeuratWrappers)
library(Banksy)

object <- RunBanksy(object,
  lambda = 0.8, verbose = TRUE,
  assay = "Spatial.008um", slot = "data", features = "variable",
  k_geom = 50
)

DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)

Idents(object) <- "banksy_cluster"
p <- SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)
p

banksy_cells <- CellsByIdentities(object)
p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
p
```

#### 6. Subset out anatomical regions
```R
cortex.coordinates <- as.data.frame(read.csv("/brahms/lis/visium_hd/final_mouse/cortex-hippocampus_coordinates.csv"))
cortex <- CreateSegmentation(cortex.coordinates)

object[["cortex"]] <- Overlay(object[["slice1.008um"]], cortex)
cortex <- subset(object, cells = Cells(object[["cortex"]]))
```

#### 7. Integration with scRNA-seq data (deconvolution)
```R
if (!requireNamespace("spacexr", quietly = TRUE)) {
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
}
library(spacexr)

# sketch the cortical subset of the Visium HD dataset
DefaultAssay(cortex) <- "Spatial.008um"
cortex <- FindVariableFeatures(cortex)
cortex <- SketchData(
  object = cortex,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(cortex) <- "sketch"
cortex <- ScaleData(cortex)
cortex <- RunPCA(cortex, assay = "sketch", reduction.name = "pca.cortex.sketch", verbose = T)
cortex <- FindNeighbors(cortex, reduction = "pca.cortex.sketch", dims = 1:50)
cortex <- RunUMAP(cortex, reduction = "pca.cortex.sketch", reduction.name = "umap.cortex.sketch", return.model = T, dims = 1:50, verbose = T)

# load in the reference scRNA-seq dataset
ref <- readRDS("/brahms/satijar/allen_scRNAseq_ref.Rds")

Idents(ref) <- "subclass_label"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$subclass_label)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- cortex[["sketch"]]$counts
cortex_cells_hd <- colnames(cortex[["sketch"]])
coords <- GetTissueCoordinates(cortex)[cortex_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
cortex <- AddMetaData(cortex, metadata = RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
cortex$first_type <- as.character(cortex$first_type)
cortex$first_type[is.na(cortex$first_type)] <- "Unknown"
cortex <- ProjectData(
  object = cortex,
  assay = "Spatial.008um",
  full.reduction = "pca.cortex",
  sketched.assay = "sketch",
  sketched.reduction = "pca.cortex.sketch",
  umap.model = "umap.cortex.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)

DefaultAssay(object) <- "Spatial.008um"

# we only ran RCTD on the cortical cells
# set labels to all other cells as "Unknown"
object[[]][, "full_first_type"] <- "Unknown"
object$full_first_type[Cells(cortex)] <- cortex$full_first_type[Cells(cortex)]
Idents(object) <- "full_first_type"

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(object)
excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
p <- SpatialDimPlot(object, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
p

plot_cell_types <- function(data, label) {
  p <- ggplot(data, aes(x = get(label), y = n, fill = full_first_type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = ifelse(n >= min_count_to_show_label, full_first_type, "")), position = position_stack(vjust = 0.5), size = 2) +
    xlab(label) +
    ylab("# of Spots") +
    ggtitle(paste0("Distribution of Cell Types across ", label)) +
    theme_minimal()
}

cell_type_banksy_counts <- object[[]] %>%
  dplyr::filter(full_first_type %in% excitatory_names) %>%
  dplyr::count(full_first_type, banksy_cluster)

min_count_to_show_label <- 20

p <- plot_cell_types(cell_type_banksy_counts, "banksy_cluster")
p

Idents(object) <- "banksy_cluster"
object$layer_id <- "Unknown"
object$layer_id[WhichCells(object, idents = c(7))] <- "Layer 2/3"
object$layer_id[WhichCells(object, idents = c(15))] <- "Layer 4"
object$layer_id[WhichCells(object, idents = c(5))] <- "Layer 5"
object$layer_id[WhichCells(object, idents = c(1))] <- "Layer 6"

# set ID to RCTD label
Idents(object) <- "full_first_type"

# Visualize distribution of 4 interneuron subtypes
inhibitory_names <- c("Sst", "Pvalb", "Vip", "Lamp5")
cell_ids <- CellsByIdentities(object, idents = inhibitory_names)
p <- SpatialDimPlot(object, cells.highlight = cell_ids, cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
p
```

#### 8. Unsupervised clustering