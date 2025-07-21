# QC2: Remove low quality cells

2025-07-22

## Remove low-quality cells

Removing low-quality cells is crucial for downstream analysis steps.
Each dataset has its own QC thresholds, and these thresholds should be determined based on reasonable evidence.
If you encounter any low-quality cells during downstream analysis, re-run the QC with adjusted thresholds!

We obtained SingleCellExperiment (SCE) objects without estimated droplets from the previous section.
In this section, we will remove low-quality cells that have large mitochondrial proportions, low RNA content, or low numbers of detected genes.

```R
library(DropletUtils)
library(dplyr)
library(scater)

set.seed(42) # for reproducibility

save_path = "QC2"
```

### Load droplet-filtered SCE objects

```R
load(file = "/BiO/data/process/QC2_data/DropletUtils_filtered_sce_1.RData")
load(file = "/BiO/data/process/QC2_data/DropletUtils_filtered_sce_2.RData")
load(file = "/BiO/data/process/QC2_data/DropletUtils_filtered_sce_3.RData")
load(file = "/BiO/data/process/QC2_data/DropletUtils_filtered_sce_4.RData")
load(file = "/BiO/data/process/QC2_data/DropletUtils_filtered_sce_5.RData")
load(file = "/BiO/data/process/QC2_data/DropletUtils_filtered_sce_6.RData")
load(file = "/BiO/data/process/QC2_data/DropletUtils_filtered_sce_7.RData")
load(file = "/BiO/data/process/QC2_data/DropletUtils_filtered_sce_8.RData")
load(file = "/BiO/data/process/QC2_data/DropletUtils_filtered_sce_9.RData")
load(file = "/BiO/data/process/QC2_data/DropletUtils_filtered_sce_12.RData")
```

### Process a single sample as an example

A SCE object contains:

- a gene-by-cell counts as a sparse matrix
- gene data (gene annotation, etc)
- cell data (sample information, experimental condition information, etc).

Gene and cell informations are stored in `rowData(SCE)` and `colData(SCE)`, respectively. Here, Ensembl gene ids are translated into gene symbol for easier further analysis.

Since there are 10 samples, we will process an example sample step-by-step and make a function `preprocess()` to apply the same process to other samples.

The example SCE object is shown below.

```R
sce <- DropletUtils_rawsce_5
sce
```

```text
class: SingleCellExperiment 
dim: 36604 4096 
metadata(1): Samples
assays(1): counts
rownames(36604): ENSG00000243485 ENSG00000237613 ... Htag2 Htag3
rowData names(3): ID Symbol Type
colnames(4096): AAACCCACATCCTAAG-1 AAACCCATCAAGCTGT-1 ...
  TTTGTTGGTGGGATTG-1 TTTGTTGTCTTCGACC-1
colData names(2): Sample Barcode
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):
```

To remove low quality cells, we use the `addPerCellQC()` function from the `scater` R package. This function calculates several quality control (QC) metrics for each cell, including:

- number of unique molecular identifiers (UMIs) per cell
- number of genes detected per cell
- percentage of UMIs assigned to mitochondrial (MT) genes

```R
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

mtgenes <- rowData(sce)[grep("^MT-", rowData(sce)$Symbol),]$Symbol
is.mito <- rownames(sce) %in% mtgenes
print(table(is.mito))
```

```text
    is.mito
    FALSE  TRUE 
    36591    13 
```

```R
sce <- addPerCellQC(
  sce,
  subsets = list(MT=mtgenes),
  percent_top = c(50, 100, 200, 500), 
  detection_limit = 5
)

sce$log10_sum = log10(sce$sum + 1)
sce$log10_detected = log10(sce$detected + 1)

sce <- sce[,sce$sum!=0]
```

We define low-quality cells with:

- &lt;500 UMIs
- &gt;15% MT gene percents
- &lt;100 detected genes

Criteria is visualized as histogram below.

```R
umi <- 500
mtpct <- 15
detect <- 100

ggplot(data.frame(colData(sce)), aes(x = log10_sum)) + 
  geom_histogram(bins = 150) +
  theme_bw() +
  geom_vline(xintercept = log10(umi), color="red", linetype="dashed")
  
ggplot(data.frame(colData(sce)), aes(x = subsets_MT_percent)) + 
  geom_histogram(bins = 150) +
  theme_bw() +
  geom_vline(xintercept = mtpct, color="red", linetype="dashed")

ggplot(data.frame(colData(sce)), aes(x = detected)) + 
  geom_histogram(bins = 150) +
  theme_bw() +
  geom_vline(xintercept = detect, color="red", linetype="dashed") + 
  ylim(0, 500)
```

![png](qc2_files/qc2_12_0.png)

![png](qc2_files/qc2_12_2.png)

![png](qc2_files/qc2_12_3.png)

Low-quality cells being filtered out can be identified by the PCA plot.

```R
filter_by_total_counts <- sce$sum > umi
filter_by_mt_percent <- sce$subsets_MT_percent < mtpct
filter_by_nfeature <- sce$detected > detect

sce <- runColDataPCA(sce, variables = list("sum", "detected", "subsets_MT_percent", "percent.top_500"))

sce$use <- (
  filter_by_total_counts &
    filter_by_mt_percent &
    filter_by_nfeature
)

plotReducedDim(sce, dimred="PCA_coldata", colour_by="use") +
  theme(plot.margin = margin(10, 10, 10, 10))

```

![png](qc2_files/qc2_14_0.png)

```R
sce <- sce[,sce$use]
sce
```

```text
class: SingleCellExperiment 
dim: 36604 1350 
metadata(1): Samples
assays(1): counts
rownames(36604): MIR1302-2HG FAM138A ... Htag2 Htag3
rowData names(3): ID Symbol Type
colnames(1350): AAACCCATCTGAACGT-1 AAACGAACAAGAGAGA-1 ...
  TTTGGTTGTCGCTCGA-1 TTTGTTGGTGGGATTG-1
colData names(15): Sample Barcode ... log10_detected use
reducedDimNames(1): PCA_coldata
mainExpName: NULL
altExpNames(0):
```

### Process other samples

Other samples can be processed in the same way, and it's baked into a function `preprocess()` below.

```R
preprocess <- function(sce){
  rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
  
  mtgenes <- rowData(sce)[grep("^MT-", rowData(sce)$Symbol),]$Symbol
  is.mito <- rownames(sce) %in% mtgenes
  
  sce <- addPerCellQC(
    sce,
    subsets = list(MT=mtgenes),
    percent_top = c(50, 100, 200, 500), 
    detection_limit = 5
  )
  
  sce$log10_sum <- log10(sce$sum + 1)
  sce$log10_detected <- log10(sce$detected + 1)
  
  sce <- sce[, sce$sum!=0]
  return(sce)
}

filtering <- function(sce,umi,mtpct,detect){
  filter_by_total_counts <- sce$sum > umi
  filter_by_mt_percent <- sce$subsets_MT_percent < mtpct
  filter_by_nfeature <- sce$detected > detect
  
  sce <- runColDataPCA(sce, variables = list("sum", "detected", "subsets_MT_percent", "percent.top_500"))
  
  sce$use <- (
    filter_by_total_counts &
      filter_by_mt_percent &
      filter_by_nfeature
  )
  
  plotReducedDim(sce, dimred="PCA_coldata", colour_by="use")
  
  sce <- sce[, sce$use]
  return(sce)
}
```

```R
umi <- 500
mtpct <- 15
detect <- 100

sce.1 <- preprocess(DropletUtils_rawsce_1)
sce.1 <- filtering(sce.1, umi, mtpct, detect)

sce.2 <- preprocess(DropletUtils_rawsce_2)
sce.2 <- filtering(sce.2, umi, mtpct, detect)

sce.3 <- preprocess(DropletUtils_rawsce_3)
sce.3 <- filtering(sce.3, umi, mtpct, detect)

sce.4 <- preprocess(DropletUtils_rawsce_4)
sce.4 <- filtering(sce.4, umi, mtpct, detect)

sce.5 <- preprocess(DropletUtils_rawsce_5)
sce.5 <- filtering(sce.5, umi, mtpct, detect)

sce.6 <- preprocess(DropletUtils_rawsce_6)
sce.6 <- filtering(sce.6, umi, mtpct, detect)

sce.7 <- preprocess(DropletUtils_rawsce_7)
sce.7 <- filtering(sce.7, umi, mtpct, detect)

sce.8 <- preprocess(DropletUtils_rawsce_8)
sce.8 <- filtering(sce.8, umi, mtpct, detect)

sce.9 <- preprocess(DropletUtils_rawsce_9)
sce.9 <- filtering(sce.9, umi, mtpct, detect)

sce.12 <- preprocess(DropletUtils_rawsce_12)
sce.12 <- filtering(sce.12, umi, mtpct, detect)
```

### Change cellbarcode name

Each cell has their own cellbarcode. To prevent errors due to the presence of the duplicate cellbarcode across samples, we will append the sample information to the cellbarcode.

Currently, cellbarcode from the first sample looks like below.

```R
head(colnames(sce.1))
```

```text
[1] "AAACCCAAGGGTGAAA-1" "AAACCCATCAGTCTTT-1" "AAACGAAGTGCGGATA-1"
[4] "AAACGAATCGTAACTG-1" "AAACGCTGTCTCACGG-1" "AAAGGATTCGTGTCAA-1"
```

```R
colnames(sce.1) <- paste0(substring(colnames(sce.1),1,16),"-L1")
colnames(sce.2) <- paste0(substring(colnames(sce.2),1,16),"-L2")
colnames(sce.3) <- paste0(substring(colnames(sce.3),1,16),"-L3")
colnames(sce.4) <- paste0(substring(colnames(sce.4),1,16),"-L4")
colnames(sce.5) <- paste0(substring(colnames(sce.5),1,16),"-L5")
colnames(sce.6) <- paste0(substring(colnames(sce.6),1,16),"-L6")
colnames(sce.7) <- paste0(substring(colnames(sce.7),1,16),"-L7")
colnames(sce.8) <- paste0(substring(colnames(sce.8),1,16),"-L8")
colnames(sce.9) <- paste0(substring(colnames(sce.9),1,16),"-L9")
colnames(sce.12) <- paste0(substring(colnames(sce.12),1,16),"-L12")
```

After appending the sample information, the cellbarcode will look like below.

```R
head(colnames(sce.1))
```

```text
[1] "AAACCCAAGGGTGAAA-L1" "AAACCCATCAGTCTTT-L1" "AAACGAAGTGCGGATA-L1"
[4] "AAACGAATCGTAACTG-L1" "AAACGCTGTCTCACGG-L1" "AAAGGATTCGTGTCAA-L1"
```

### Adding sample information as metadata

For downstream analysis, we will combine all SCE objects add metadata about their conditions (`CONTROL`, `ASYMPTOMATIC` or `SYMPTOMATIC`)

```R
sce <- cbind(sce.1, sce.2, sce.3, sce.4, sce.5, sce.6, sce.7, sce.8, sce.9, sce.12)
sce$library <- "NA"
sce$library[grep("20094_0001_A_B", sce$Sample)] <- "L1"
sce$library[grep("20094_0002_A_B", sce$Sample)] <- "L2"
sce$library[grep("20094_0003_A_B", sce$Sample)] <- "L3"
sce$library[grep("20094_0004_A_B", sce$Sample)] <- "L4"
sce$library[grep("20094_0005_A_B", sce$Sample)] <- "L5"
sce$library[grep("20094_0006_A_B", sce$Sample)] <- "L6"
sce$library[grep("20094_0007_A_B", sce$Sample)] <- "L7"
sce$library[grep("20094_0008_A_B", sce$Sample)] <- "L8"
sce$library[grep("20094_0009_A_B", sce$Sample)] <- "L9"
sce$library[grep("20094_0012_A_B", sce$Sample)] <- "L12"

table(sce$library)
```

```text
      L1  L12   L2   L3   L4   L5   L6   L7   L8   L9 
    1823 3293 2166 1835  833 1350  324  419 1345  678 
```

```R
sce$Condition <- "SYMPTOMATIC"
sce$Condition[sce$library == "L12"] <- "CONTROL"
sce$Condition[sce$library %in% c("L4","L5","L6")] <- "ASYMPTOMATIC"

table(sce$Condition)
```

```text
    ASYMPTOMATIC      CONTROL  SYMPTOMATIC 
            2507         3293         8266 
```

## References

Butler, A., Hoffman, P., Smibert, P., Papalexi, E. & Satija, R. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat. Biotechnol. 36, 411-420 (2018).

McCarthy, D. J., Campbell, K. R., Lun, A. T. & Wills, Q. F. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. Bioinformatics 33, 1179-1186 (2017).

Pekayvaz, K., Leunig, A., Kaiser, R. et al. Protective immune trajectories in early viral containment of non-pneumonic SARS-CoV-2 infection. Nat Commun 13, 1018 (2022).
