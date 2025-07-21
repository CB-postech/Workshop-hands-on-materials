QC1-DropletUtils
================
2025-07-18

## QC thresholds should be vary data to data
## These various threshold can be determined by visual inspection (drawing histogram, ...)

### QC1 Remove Empty Droplets
Many droplets are empty (mainly consist of ambient RNA) <br>
We have to filter out these droplets <br>
1st) Construct an ambient pool by summing up UMIs from low-count droplet <br>
2nd) Test the possibility that ambient RNA accidentally generated this droplet profile <br>
These process are wrapped up by package DropletUtils


```R
library(Seurat)
library(scater)
library(DropletUtils)
library(scran)
library(ggplot2)

set.seed(42) # for reproducibility

rdatadir = '/BiO/data/QC/'
save_path = 'QC1'
```

### load


```R
rawsce_1 <- read10xCounts(paste0(rdatadir, "20094_0001_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_2 <- read10xCounts(paste0(rdatadir, "20094_0002_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_3 <- read10xCounts(paste0(rdatadir, "20094_0003_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_4 <- read10xCounts(paste0(rdatadir, "20094_0004_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_5 <- read10xCounts(paste0(rdatadir, "20094_0005_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_6 <- read10xCounts(paste0(rdatadir, "20094_0006_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_7 <- read10xCounts(paste0(rdatadir, "20094_0007_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_8 <- read10xCounts(paste0(rdatadir, "20094_0008_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_9 <- read10xCounts(paste0(rdatadir, "20094_0009_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_12 <- read10xCounts(paste0(rdatadir, "20094_0012_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
```

### check


```R
rawsce_1
```


```R
# class: SingleCellExperiment 
# dim: 36604 6794880 
# metadata(1): Samples
# assays(1): counts
# rownames(36604): ENSG00000243485 ENSG00000237613 ... Htag2 Htag3
# rowData names(3): ID Symbol Type
# colnames: NULL
# colData names(2): Sample Barcode
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
```

### Run DropletUtils


```R
rawsce_list <- c('rawsce_1', 'rawsce_2', 'rawsce_3', 'rawsce_4', 'rawsce_5',
                 'rawsce_6', 'rawsce_7', 'rawsce_8', 'rawsce_9', 'rawsce_12')
```


```R

for (i in 1:length(rawsce_list)) {
  rawsce <- get(rawsce_list[i])
  br.out <- barcodeRanks(counts(rawsce))
  
  e.out <- emptyDrops(counts(rawsce))  ## Cells that have UMI counts lower than 100 (by defualt) are empty cells.
  print(table(Sig=e.out$FDR <= 0.05, Limited=e.out$Limited))

  is.cell <- e.out$FDR <= 0.05
  print(sum(is.cell, na.rm=TRUE))

  p <- ggplot(data.frame(br.out), aes(x = rank, y = total)) + 
    geom_point() + 
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    ggtitle(rawsce_list[i]) +
    theme_classic()
	p <- p + 
		geom_hline(aes(yintercept = min(br.out$fitted, na.rm = TRUE),  color = "FDR_0.05"), linetype = "dashed") +
		geom_hline(aes(yintercept = as.numeric(metadata(br.out)$knee), color = "knee"), linetype = "dashed") +
		geom_hline(aes(yintercept = as.numeric(metadata(br.out)$inflection), color = "inflection"), linetype = "dashed") +
		scale_color_manual(name = "Cutoffs",
											values = c("FDR_0.05" = "red",
																	"knee" = "dodgerblue",
																	"inflection"= "forestgreen")) +
		theme(legend.position = "bottom")
  ggsave(filename = paste0(save_path, '/DropletUtils_', rawsce_list[i], '.png'), plot = p, width = 6, height = 6)
  
  colnames(rawsce) = colData(rawsce)$Barcode
  rawsce <- rawsce[,which(e.out$FDR <= 0.05)]
  
  assign(paste0('DropletUtils_', rawsce_list[i]), rawsce)
}
```
<img width="1800" height="1800" alt="image" src="https://github.com/user-attachments/assets/25285f61-2a6a-4c6f-99ca-636bfdb49a7e" />


```R
##        Limited
## Sig     FALSE  TRUE
##   FALSE 65265     0
##   TRUE   5709  6541
## [1] 12250

##        Limited
## Sig     FALSE  TRUE
##   FALSE 75151     0
##   TRUE   1406  5899
## [1] 7305

##        Limited
## Sig     FALSE  TRUE
##   FALSE 68613     0
##   TRUE   660   5843
## [1] 6503

##        Limited
## Sig     FALSE  TRUE
##   FALSE 30934     0
##   TRUE    561  3931
## [1] 4492

##        Limited
## Sig     FALSE  TRUE
##   FALSE 10689     0
##   TRUE    428  3671
## [1] 4099

##        Limited
## Sig     FALSE  TRUE
##   FALSE 25095     0
##   TRUE    427  2209
## [1] 2636

##        Limited
## Sig     FALSE  TRUE
##   FALSE 14106     0
##   TRUE    196  1806
## [1] 2002

##        Limited
## Sig     FALSE  TRUE
##   FALSE 68921     0
##   TRUE    869  4266
## [1] 5135

##        Limited
## Sig     FALSE  TRUE
##   FALSE 67020     0
##   TRUE    549  2325
## [1] 2874

##        Limited
## Sig     FALSE  TRUE
##   FALSE 60514     0
##   TRUE    689  6462
## [1] 7151
```

### save


```R

save(DropletUtils_rawsce_1, file = paste0(save_path, '/DropletUtils_filtered_sce_1.RData'))
save(DropletUtils_rawsce_2, file = paste0(save_path, '/DropletUtils_filtered_sce_2.RData'))
save(DropletUtils_rawsce_3, file = paste0(save_path, '/DropletUtils_filtered_sce_3.RData'))
save(DropletUtils_rawsce_4, file = paste0(save_path, '/DropletUtils_filtered_sce_4.RData'))
save(DropletUtils_rawsce_5, file = paste0(save_path, '/DropletUtils_filtered_sce_5.RData'))
save(DropletUtils_rawsce_6, file = paste0(save_path, '/DropletUtils_filtered_sce_6.RData'))
save(DropletUtils_rawsce_7, file = paste0(save_path, '/DropletUtils_filtered_sce_7.RData'))
save(DropletUtils_rawsce_8, file = paste0(save_path, '/DropletUtils_filtered_sce_8.RData'))
save(DropletUtils_rawsce_9, file = paste0(save_path, '/DropletUtils_filtered_sce_9.RData'))
save(DropletUtils_rawsce_12, file = paste0(save_path, '/DropletUtils_filtered_sce_12.RData'))
```

### Reference
Lun, A. T., Riesenfeld, S., Andrews, T., Gomes, T., & Marioni, J. C. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome biology, 20(1), 1-9.

Pekayvaz, K., Leunig, A., Kaiser, R., Joppich, M., Brambs, S., Janjic, A., … & Nicolai, L. (2022). Protective immune trajectories in early viral containment of non-pneumonic SARS-CoV-2 infection. Nature communications, 13(1), 1-21.

