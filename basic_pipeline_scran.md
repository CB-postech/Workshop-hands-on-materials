# Basic Pipeline
---------------------------------
2025.07.18
1. Normalization
2. Batch-aware HVG selection
3. Dimensional Reduction
4. Batch-Correciton
5. Annotation
6. celltypist prediction
7. DEG Analysis


```R
library(scran)
library(scater)
library(ggplot2)
set.seed(42) # for reproducibility

save_path = "basic_pipeline_scran"
setwd(save_path)
```


### 0. load data
data from https://www.nature.com/articles/s41591-023-02327-2<br>
You can access the data via CellxGene (https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293)<br>
In this example, we will use the sampled data


```R
count_matrix <- read.csv("/BiO/data/process/basic_pipeline_data/HLCA_pulmonary_fibrosis_immune_raw.csv", row.names = 1)
meta.data <- read.csv("/BiO/data/process/basic_pipeline_data/HLCA_pulmonary_fibrosis_immune_meta.csv", row.names = 1)

sce <- SingleCellExperiment(list(counts=count_matrix))
```


### 1. Normalization


```R
# scran style
# To remove cell-specific biases, cells are clustered using quickCluster() and cell-specific size factors are calculated using computeSumFactors() of scran R package. 
# Raw counts of each cell are divided by cell-specific size factor and log2-transformed with a pseudocount of 1.

clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
sce.norm <- logNormCounts(sce,  pseudo_count = 1)

library(Seurat)
so <- as.Seurat(sce.norm,
                    counts = "counts",
                    data = "logcounts")

so$disease = meta.data[Cells(so), 'disease']
so$study = meta.data[Cells(so), 'study']

# change assay name
so[['RNA']] = so[['originalexp']]
DefaultAssay(so) = 'RNA'
so[['originalexp']] = NULL
```

### 2. Batch-aware Feature selection
In this data, there are multiple studies. <br>
Therefore, we will select HVG with consideration of the source of data (study). <br>
batch-aware HVG selection is implemented in Scanpy python library <br>


```R
table(so$study)
```


    
     Banovich_Kropski_2020          Kaminski_2020 Misharin_Budinger_2018 
                      2338                   2328                    810 
             Sheppard_2020 
                      2306 


```R
source('/BiO/data/batch_aware_in_seurat.R')
```


```R
batch_key = 'study'
nHVG = 2000
batch_aware_in_seurat(so, batch_key = batch_key, nHVG = nHVG, conda_env = "/BiO/prog/miniforge3/envs/QC", save_path = save_path)
```


```R
HVG = read.csv(paste0(save_path, "/hvg_", nHVG, "_", batch_key, ".csv"))
# > head(HVG)
#         gene
# 1        A2M
# 2       ABAT
# 3      ABCA1
# 4      ABCA2
# 5      ABCA6
# 6      ABCB1
```

### 3. Dimensionality Reduction
Normalized data is scaled and principal components (PCs) are calculated by a gene-by-cell matrix with HVGs.


```R
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)

so <- RunPCA(so, features = HVG$gene)

p <- ElbowPlot(so, ndims = 50)
ggsave(p, filename = paste0(save_path, "/elbow_plot.png"), width = 4, height = 4)

PCs <- 8
so <- FindNeighbors(so, dims = 1:PCs)
so <- FindClusters(so, resolution = 0.5)

so <- RunUMAP(so, dims = 1:PCs)
# so <- RunTSNE(so, dims = 1:PCs)

DimPlot(so, group.by = 'study')

# save as png
p <- DimPlot(so, group.by = 'study')
ggsave(p, filename = paste0(save_path, "/umap_study.png"), width = 8, height = 6)
```
<img width="840" height="840" alt="output_6_3" src="https://github.com/user-attachments/assets/5289c2ae-ccfd-4e43-b106-1c6320650061" />
<img width="840" height="840" alt="output_7_0" src="https://github.com/user-attachments/assets/b57d9f63-2b38-4c3b-a6ea-277ab0ca3ad3" />


### 4. Batch Correction by Harmony
ref) Korsunsky, I., Millard, N., Fan, J. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289–1296 (2019). https://doi.org/10.1038/s41592-019-0619-0 <br>
https://portals.broadinstitute.org/harmony/articles/quickstart.html


```R
library(harmony)
so <- RunHarmony(so, 'study')
so <- FindNeighbors(so, reduction = "harmony")
so <- FindClusters(so, resolution = 0.5) 
so <- RunUMAP(so, dims = 1:PCs, reduction = 'harmony', reduction.name = 'umap.harmony') # use same dimension number as before
# so <- RunTSNE(so, dims = 1:PCs, reduction = 'harmony', reduction.name = 'tsne.harmony')

# save as png
p <- DimPlot(so, group.by = 'study', reduction = 'umap.harmony')
ggsave(p, filename = paste0(save_path, "/umap_harmony_study.png"), width = 8, height = 6)
p <- DimPlot(so, group.by = 'seurat_clusters', reduction = 'umap.harmony', label = TRUE)
ggsave(p, filename = paste0(save_path, "/umap_harmony_seurat_clusters.png"), width = 8, height = 6)
```

```R
DimPlot(so, group.by = 'study', reduction = 'umap.harmony')
DimPlot(so, group.by = 'seurat_clusters', reduction = 'umap.harmony', label = TRUE)
```
<img width="840" height="840" alt="output_9_1" src="https://github.com/user-attachments/assets/fb1c8b82-dcbf-49ee-92e0-f2c8dd4bfd08" />
<img width="840" height="840" alt="output_9_0" src="https://github.com/user-attachments/assets/3853777d-c0cf-4b5b-a322-108ecbfb4b35" />

### 5. Celltypist Prediction


```R
source("/BiO/data/celltypist_in_seurat.R")
celltypist_in_seurat(so, conda_env = "/BiO/prog/miniforge3/envs/QC", save_path = save_path, model_path = '/BiO/data/Immune_All_High.pkl')
```
<img width="461" height="255" alt="CellTypist_dotplot_celltypist_dotplot_majority_voting" src="https://github.com/user-attachments/assets/50cad2ce-4adc-4ba2-b885-5496f3e11d20" />
<img width="552" height="759" alt="CellTypist_dotplot_celltypist_dotplot_predicted_labels" src="https://github.com/user-attachments/assets/4fa6258b-0069-4ede-a4ee-9ebe9959a1e1" />


### 6. Marker gene expression visualization


```R
marker.genes<- list(T.cell = c('CD3D', 'CD3E'),
                    NK.cell = c('KLRD1', 'NKG7', 'GNLY'),
                    Macrophage = c('FABP4', 'APOE', 'MARCO'),
                    Monocyte = c('S100A12', 'FCN1'))
p <- FeaturePlot(so, features = as.vector(unlist(marker.genes)), ncol = 4, reduction = 'umap.harmony')
ggsave(p, filename = paste0(save_path, "/feature_plot_markers.png"), width = 13, height = 9)

p <- DotPlot(so, features = unlist(marker.genes), group.by = 'seurat_clusters')
p <- p + theme(axis.text.x = element_text(angle = 90))
ggsave(p, filename = paste0(save_path, "/dot_plot_markers.png"), width = 10, height = 6)

so$celltype = as.character(so$seurat_clusters)
so$celltype[so$celltype %in% c(0, 5)] <- 'T cell'
so$celltype[so$celltype %in% c(3)] <- 'NK cell'
so$celltype[so$celltype %in% c(1, 4, 7)] <- 'Macrophage'
so$celltype[so$celltype %in% c(2, 6)] <- 'Monocyte'
so$celltype = factor(so$celltype, levels = c('T cell', 'NK cell', 'Macrophage', 'Monocyte'))

p <- DimPlot(so, group.by = 'celltype', reduction = 'umap.harmony', label = TRUE)
ggsave(p, filename = paste0(save_path, "/umap_harmony_celltype.png"), width = 8, height = 6)

p <- DotPlot(so, features = as.vector(unlist(marker.genes)), group.by = 'celltype')
p <- p + theme(axis.text.x = element_text(angle = 90))
ggsave(p, filename = paste0(save_path, "/dot_plot_celltype_markers.png"), width = 10, height = 6)

# saveRDS(so, file = paste0(save_path, "/HLCA_pulmonary_fibrosis_immune.rds"))
```



```R
FeaturePlot(so, features = as.vector(unlist(marker.genes)), ncol = 4, reduction = 'umap.harmony')
DotPlot(so, features = unlist(marker.genes), group.by = 'seurat_clusters') + theme(axis.text.x = element_text(angle = 90))
DimPlot(so, group.by = 'celltype', reduction = 'umap.harmony', label = TRUE)
DotPlot(so, features = as.vector(unlist(marker.genes)), group.by = 'celltype') + theme(axis.text.x = element_text(angle = 90))
```
<img width="840" height="840" alt="output_12_4" src="https://github.com/user-attachments/assets/e2ec058d-d006-4b5c-876f-29d6e0565f0d" />
<img width="840" height="840" alt="output_12_3" src="https://github.com/user-attachments/assets/99d4b8fb-b9f0-4c84-9afa-dc92750c2380" />
<img width="840" height="840" alt="output_12_1" src="https://github.com/user-attachments/assets/a7ca2a43-385b-43c6-8c3a-baebad3e12ef" />
<img width="840" height="840" alt="output_12_0" src="https://github.com/user-attachments/assets/1fb8131b-9eeb-42fb-b9e5-8504c4d908fb" />


### 7. Patient-aware DEG selection


```R
donor_id <- read.csv('/BiO/data/HLCA_pulmonary_fibrosis_donor_subset.csv', row.names = 1)
so$donor.id <- donor_id[Cells(so), ]
```


```R
library(EnhancedVolcano)

Tcells.control = Cells(so)[so$celltype == 'T cell' & so$disease == 'normal']
Tcells.pulmonary = Cells(so)[so$celltype == 'T cell' & so$disease == 'pulmonary fibrosis']

DEG.vanilla <- FindMarkers(so, ident.1 = Tcells.pulmonary, ident.2 = Tcells.control)
DEG.patient.correction <- FindMarkers(so, ident.1 = Tcells.pulmonary, ident.2 = Tcells.control, test.use = 'MAST', latent.vars = 'donor.id')

p <- EnhancedVolcano(DEG.vanilla,
    lab = rownames(DEG.vanilla),
    x = 'avg_log2FC',
    y = 'p_val_adj', 
    FCcutoff = log2(1.5),
    pCutoff = 0.05)
ggsave(p, filename = paste0(save_path, "/volcano_Tcells_vanilla.png"), width = 8, height = 6)

p <- EnhancedVolcano(DEG.patient.correction,
    lab = rownames(DEG.patient.correction),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    FCcutoff = log2(1.5),
    pCutoff = 0.05)
ggsave(p, filename = paste0(save_path, "/volcano_Tcells_patient_correction.png"), width = 8, height = 6)
```

<img width="840" height="840" alt="output_15_1" src="https://github.com/user-attachments/assets/194434cf-e68a-40a4-baee-6484a758b8c9" />
<img width="840" height="840" alt="output_15_0" src="https://github.com/user-attachments/assets/ab1b6a9c-481a-4883-8b41-0988e0a9b8c5" />
