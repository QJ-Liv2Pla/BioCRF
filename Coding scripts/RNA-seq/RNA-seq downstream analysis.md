```shell
source /opt/annaconda/annaconda3/bin/activate
conda activate env_nf
```

## common features of RNA-seq count data:
* a low number of counts associated with a large proportion of genes
* a long right tail due to the lack of any upper limit for expression
* large dynamic range.

## Modeling count data
Count data in general can be modeled with various distributions:
* Binomial distribution: Gives you the probability of getting a number of heads upon tossing a coin a number of times. Based on discrete events and used in situations when you have a certain number of cases.

* Poisson distribution: For use, when the number of cases is very large (i.e. people who buy lottery tickets), but the probability of an event is very small (probability of winning). The Poisson is similar to the binomial, but is based on continuous events. Appropriate for data where mean == variance.

* Negative binomial distribution: An approximation of the Poisson, but has an additional parameter that adjusts the variance independently from the mean.



### Load R packages
```r
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
```
### Load data
```shell
ls -d */ | grep 'RNA' | xargs -I {} cp -r {} /data1/Analysis/NGS_analysis_20240820/RNA_seq_downstream_analysis
cp /hdd/202408_NGS_raw/Aligned_data_20240826/RNAseq_analysis/star_salmon/tx2gene.tsv /data1/Analysis/NGS_analysis_20240820/RNA_seq_downstream_analysis
```

```r
## List all directories containing data  
samples <- list.files(path = "./data", full.names = T, pattern="salmon$")

## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")

## Since all quant files have the same name it is useful to have names for each element
names(files) <- str_replace(samples, "./data/", "") %>% 
                str_replace(".salmon", "")
```


```r
# Load the annotation table for GrCh38
tx2gene <- read.delim("")
tx2gene %>% View()
```

```r
?tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "ensgene")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE)
attributes(txi)
```  

```r
txi$counts %>% View()
# Write the counts to an object
data <- txi$counts %>% 
  round() %>% 
  data.frame()
```

```r
## Create a sampletable/metadata
sampletype <- factor(rep("temp37",3),rep("temp39",3),rep("temp42",3))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))
```

### Check RNA-seq count distribution
```r
ggplot(data) +
  geom_histogram(aes(x = temp37_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
```

```r
mean_counts <- apply(data[,6:8], 1, mean) #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(data[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
        geom_point(aes(x=mean_counts, y=variance_counts)) + 
        scale_y_log10(limits = c(1,1e9)) +
        scale_x_log10(limits = c(1,1e9)) +
        geom_abline(intercept = 0, slope = 1, color="red")
```

```r
### Check that sample names match in both files
all(colnames(txi$counts) %in% rownames(meta))
all(colnames(txi$counts) == rownames(meta))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)
```