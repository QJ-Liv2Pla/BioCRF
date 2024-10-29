# Task
> Dinopore for A to I editing identification in ONT sequencing data
#
# Step 1. Environment configuration
## (Option)Use [Dinopore Docker file](https://codeocean.com/capsule/4038948/tree/v1) to set up environment directly

## Otherwise set up the environment step by step
## 1. Conda install
> ```Annaconda``` provide a clean and isolate environment for different tasks. Packages/Softwares can be easily installed through [```conda-forge```](https://anaconda.org/conda-forge). This is one of the best practice to avoid environment conflicts.

Follow the instruction on [conda install website](https://docs.anaconda.com/miniconda/), miniconda works fine.
```shell
# Download
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
# Install
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
# Activate
source ~/miniconda3/bin/activate
```

## 2. REQUISITE SOFTWARE
### 2.1 List
Guppy_basecaller 3.2.4 (Please modify S1.Basecall_map_nanopolish.sh to point to it)

Third party software included
- Graphmap2		 			0.6.3 
- Sam2tsv (part of Jvarkit)	34d8e7f7 
- Picard					2.21.6

conda env
- python		3.8.5
- h5py			2.10.0
- nanopolish	0.11.1
- pillow		8.3.1
- pyyaml		5.4.1
- requests		2.26.0
- samtools		1.9
- scipy			1.7.1

R Packages (R version 4.1)
- Matrix		1.3.4
- R.utils		2.11.0
- Rcpp			1.0.7
- abind			1.4.5
- caret			6.0.89
- data.table	1.14.2
- doParallel	1.0.16
- ff			4.0.4
- foreach		1.5.1
- keras			>= 2.3.0
- multiROC		1.1.1
- optparse		1.6.6
- pacman		0.5.1
- plyr			1.8.6
- pracma		2.3.3
- scales		1.1.1
- tensorflow	>= 2.3.0
- tidyverse		1.3.1
- usefun		0.4.8
- zoo			1.8.9

### 2.2 INSTALLATION
>To install them one by one
#### 2.2.1 Conda env
```shell
# Create a new environment
conda create -n Dinopore python=3.8.5 h5py=2.10.0 nanopolish=0.11.1 pillow=8.3.1 pyyaml=5.4.1 requests=2.26.0 samtools=1.9 scipy=1.7.1 r-base=4.1.5 
# Activate the new env
conda activate Dinopore
# Install following packages through conda-forge
conda install conda-forge::r-nloptr --yes
conda install conda-forge::r-caret --yes
conda install conda-forge::r-data.table --yes
conda install conda-forge::r-doparallel --yes
conda install conda-forge::r-tidyverse --yes
conda install conda-forge::openjdk=17.0.1 --yes

# Enter R
R
```
#### 2.2.2 R
```R
# Use R command
install.packages("tensorflow")
library(tensorflow)
install_tensorflow(envname = "Dinopore", version = "2.3.0")
install.packages(c('Matrix','R.utils','Rcpp','abind','ff','foreach','keras','multiROC','optparse','pacman','plyr','pracma','scales','usefun','zoo'))
quit()
n
```
#### 2.2.3 Software
[Guppy_basecaller installation guide](https://github.com/kishwarshafin/lab-notes/blob/master/basecalling_with_guppy.md) (Please modify S1.Basecall_map_nanopolish.sh to point to it)

[Graphmap2 installation guide](https://github.com/lbcb-sci/graphmap2)

[Picard installation guide](https://github.com/broadinstitute/picard)


# Step2 Run Dinopore

To convert pod5 to fast5, create a new conda env to install pod5
```shell
conda create -n pod5 python=3.10
conda activate pod5
conda install h5py=3.11.0 --yes
pip install pod5
pod5 convert to_fast5 --help
pod5 convert to_fast5 -r /data3/Nanopore_analysis/Dinopore/input_data/pod5 --output /data3/Nanopore_analysis/Dinopore/input_data/fast5
ls /data3/Nanopore_analysis/Dinopore/input_data/fast5
```

```shell
# Open a terminal
screen -S Dinopore
```
#### Use dorado to speed up alignment Or use Dinopore script
```shell
#! /bin/bash

#SBATCH -J dorado_sup
#SBATCH -p gpu1
#SBATCH -N 1
#SBATCH --gres=gpu:3
#SBATCH --cpus-per-gpu=20

/cluster/home/hongyanhong/soft/dorado-0.7.2-linux-x64/bin/dorado basecaller  sup /cluster/home/hongyanhong/projects/CDKN1A_RAB44/enrich_sample_forall_240725/2024-5-18-YHY/2024_5_18_YHY/no_sample/20240518_2127_MN28697_FAZ02462_430ad21e/pod5  --kit-name SQK-PCB111-24 -v \
 > fastq_dorado_sup/calls.bam 
```

#### 
> Before use, edit every .sh and .R script to modify file/script path. e.g. ```./code/``` instead of ```/code/``` after you git clone the Dinopore script.
```shell
flowcell=FLO-PRO004RA
kit=SQK-RNA004
# Step 1 - Basecall fast5 -> map to genome reference -> run nanopolish to extract signal

exptdir=/data3/Nanopore_analysis/Dinopore/input_data/
ref=/data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
numcore=16
delfastq='y'

# Step 2 - Process data: convert bam file to tsv and combine nanopolish into single signal for each 5-mer of a read

exptdir=
ref=
numcore=

# Step 3 - Combine bam-tsv file and combined nanopolish file to generate raw features table

exptdir=
numcore=
agggrp='PAW01419'

# Step 4 - Aggregate features of reads into positions
exptdir=
numcore=
agggrp=

# Step 5 - Transform 1D into 2D data + Label data (class 0, 1 and 2 for unmodified, Inosine and SNP AG)
exptdir=
numcore=
agggrp=
classref=/data/xen_s9_r1_50k/groundtruth_classification.txt # Use proper groundtruth data

# Step 6
exptdir=
numcore=
agggrp=
classref=
```
#### Run Dinopore
```Shell
# e.g.
cd /data3/Nanopore_analysis/Dinopore
bash mainscript.sh -e /data3/Nanopore_analysis/Dinopore/input_data/ -r /data1/Data/Reference/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -n 16 -g 'PAW01419' -c groundtruth_classification.txt
```