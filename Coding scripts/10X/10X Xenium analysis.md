# Xenium Ranger
## Xenium Ranger is a set of analysis pipelines that process Xenium In Situ Gene Expression data. It includes the following pipelines:
* Relabel transcripts with ```xeniumranger relabel```
* Resegment cells with latest 10x segmentation algorithms with ```xeniumranger resegment```
* Import your own segmentation data to assign transcripts to cells with ```xeniumranger import-segmentation```
* Rename region or cassette names with ```xeniumranger rename```

## Download & install Xeniumranger
```shell
cd /opt
sudo mkdir Xeniumranger
cd Xeniumranger
curl -o xeniumranger-2.0.1.tar.gz "https://cf.10xgenomics.com/releases/xeniumranger/xeniumranger-2.0.1.tar.gz?Expires=1725300474&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=BITYA48Bx-M4Md6K-RFIuAPhKdeZbpNgq6izDUNPLCO0OrPT5wY9E6tK8w4xK7nX2ER2YrN0faTrb13IjNztjZ6uYHfa~dfAZTlC2S361tWVxhG9l192jn9cniW7YhNPhVTCXl2Evkpu5U~4pVTV4O4SaKwYRwjNYnPmXCG6fMvC-Y6fbzQEsU3JrDGloxoYwU22vtD789yCnfqlCGaHU128-ypJ5loZDwa7wMalKXFXmuPMP9Q8vOCtBVof~MubPFaHN6ESy59WNCckGyxFm1ABRqTYQdYnjaZG-EhFqSLjtIZdMNPsLQ4X4yZP0QYPYF~SDIc5ciNeN2nAASijsA__"
tar -xzvf xeniumranger-2.0.1.tar.gz
export PATH=/opt/Xeniumranger/xeniumranger-xenium2.0.1:$PATH
```

## (Optional) Relabel Decoded Transcripts
The ```relabel``` pipeline allows you to change the gene labels applied to decoded transcripts. If the wrong gene_panel.json file was selected on the Xenium Analyzer, you can use this pipeline to correct the gene labels after the instrument run completes. 

```shell
xeniumranger relabel --id=relabel-demo \
                     --xenium-bundle=/path/to/xenium/files \
                     --panel=/opt/xeniumranger-xenium2.0/lib/json/definitions/panel_designer/panels/xenium_human_brain_gene_expression_panel.json \
                     --localcores=32 \
                     --localmem=128
```

## (Optional) Relabel Decoded Transcripts
The ```rename``` pipeline allows you to change the sample ```region_name``` and ```cassette_name``` throughout all the Xenium Onboard Analysis output files that contain this information.
```shell
xeniumranger rename --id=rename-demo \
                     --xenium-bundle=/path/to/xenium/files \
                     --region-name="replicate 1" \
                     --cassette-name="cassette 1"
```
## (Optional) Rerun Xenium Onboard Analysis Segmentation
The ```resegment``` pipeline allows you to generate a new segmentation of the morphology image space
```shell
xeniumranger resegment --id=demo \
                       --xenium-bundle=/path/to/xenium/files \
                       --localcores=32 \
                       --localmem=128
```

## Import Nucleus and Cell Segmentation Results
The ```import-segmentation``` pipeline allows you to specify 2D nuclei and/or cell segmentation results to use for assigning transcripts to cells and recalculate all Xenium Onboard Analysis (XOA) outputs that depend on segmentation.
```shell
# Only change nuclear expansion
xeniumranger import-segmentation --id=demo \
                                 --xenium-bundle=/path/to/xenium/files \
                                 --nuclei=/path/to/cells.zarr.zip \
                                 --localcores=32 \
                                 --localmem=128
                    
# Nucleus-only count matrix
xeniumranger import-segmentation --id=demo \
                                 --xenium-bundle=/path/to/xenium/files \
                                 --nuclei=/path/to/cells.zarr.zip \
                                 --expansion-distance=0 \
                                 --localcores=32 \
                                 --localmem=128

# Polygon input
xeniumranger import-segmentation --id=polygon-demo \
                                 --xenium-bundle=/path/to/xenium/files \
                                 --nuclei=polygon.geojson \
                                 --cells=polygon.geojson \
                                 --units=pixels \
                                 --localcores=32 \
                                 --localmem=128
                    
# Mask input
xeniumranger import-segmentation --id=mask-demo-nuclei \
                                 --xenium-bundle=/path/to/xenium/files \
                                 --nuclei=nuclei.tif \
                                 --cells=cells.npy \
                                 --localcores=32 \
                                 --localmem=128

# Transcript assignment input
xeniumranger import-segmentation --id=baysor-demo \
                                 --xenium-bundle=/path/to/xenium/files \
                                 --transcript-assignment=segmentation.csv \
                                 --viz-polygons=segmentation_polygons.json \
                                 --units=microns \
                                 --localcores=32 \
                                 --localmem=128

# Transformation matrix input
xeniumranger import-segmentation --id=image-transform-demo \
                                 --xenium-bundle=/path/to/xenium/files \
                                 --coordinate-transform=demo_imagealignment.csv \
                                 --nuclei=qupath.geojson \
                                 --cells=qupath.geojson \
                                 --units=pixels \
                                 --localcores=32 \
                                 --localmem=128
```
## Output
* Analysis summary: check metrics and plots on the Cell Segmentation tab. There are new metrics for imported segmentation results (present in ```analysis_summary.html``` and ```metrics_summary.csv```).
* All secondary analysis files
* All cell, nuclei, and transcript files
