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
library(Seurat)
library(harmony)
library(ggplot2)
set.seed(42) # for reproducibility

save_path = "basic_pipeline"
```


### 0. load data
data from https://www.nature.com/articles/s41591-023-02327-2<br>
You can access the data via CellxGene (https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293)<br>
In this example, we will use the sampled data


```R
count_matrix <- read.csv("/BiO/data/process/basic_pipeline_data/HLCA_pulmonary_fibrosis_immune_raw.csv", row.names = 1)
meta.data <- read.csv("/BiO/data/process/basic_pipeline_data/HLCA_pulmonary_fibrosis_immune_meta.csv", row.names = 1)

so <- CreateSeuratObject(counts = count_matrix, meta.data = meta.data, assay = 'RNA', min.cells = 0, min.features = 0, project = 'HLCA_Pulmonary_Fibrosis_immune')
# genes are in rows, cells are in columns

# so stand for 's'eurat 'o'bject
# In this example, we use filtered data, so set min.cells and min.features to 0 (no filtering)

# > head(so, n = 3)
#                                                      orig.ident nCount_RNA nFeature_RNA            disease                 study
# F01173_GCTGGGTTCCTGTAGA_haberman HLCA_Pulmonary_Fibrosis_immune       5525         1877 pulmonary fibrosis Banovich_Kropski_2020
# F00431_CTAGAGTCATGCCACG_haberman HLCA_Pulmonary_Fibrosis_immune       2784         1017 pulmonary fibrosis Banovich_Kropski_2020
# F01172_AGTAGTCGTCCGACGT_haberman HLCA_Pulmonary_Fibrosis_immune       1617         1012 pulmonary fibrosis Banovich_Kropski_2020
```

    Warning message:
    “Data is of class data.frame. Coercing to dgCMatrix.”


### 1. Normalization


```R
# seurat style
so <- NormalizeData(so)
```

    Normalizing layer: counts
    

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
# 2       AACS
# 3 AARS.AARS1
# 4      ABCA1
# 5     ABCA13
# 6      ABCA2
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
<img width="840" height="840" alt="output_6_0" src="https://github.com/user-attachments/assets/5194a94a-739d-489f-83ef-b8f73fc59e6a" />
<img width="840" height="840" alt="output_5_3" src="https://github.com/user-attachments/assets/8e7111d0-dcd0-4672-b9ff-67ce22cabdbe" />


### 4. Batch Correction by Harmony
ref) Korsunsky, I., Millard, N., Fan, J. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289–1296 (2019). https://doi.org/10.1038/s41592-019-0619-0 <br>
https://portals.broadinstitute.org/harmony/articles/quickstart.html


```R
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
<img width="840" height="840" alt="output_7_4" src="https://github.com/user-attachments/assets/71648fa2-4c2c-482a-90be-128312bcb793" />
<img width="840" height="840" alt="output_7_3" src="https://github.com/user-attachments/assets/8946998f-6af5-4d7d-9914-e746a6a5f594" />

### 5. Celltypist Prediction


```R
source("/BiO/data/celltypist_in_seurat.R")
celltypist_in_seurat(so, conda_env = "/BiO/prog/miniforge3/envs/QC", save_path = save_path, model_path = '/BiO/data/Immune_All_High.pkl')
```
<img width="490" height="255" alt="CellTypist_dotplot_celltypist_dotplot_majority_voting (1)" src="https://github.com/user-attachments/assets/29e4fd8a-2ba5-476a-96b0-311d552f86d4" />
<img width="461" height="255" alt="CellTypist_dotplot_celltypist_dotplot_majority_voting" src="https://github.com/user-attachments/assets/4eb6995f-ac5e-4b66-b482-5d861ed0b036" />



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
so$celltype[so$celltype %in% c(1, 6, 8)] <- 'T cell'
so$celltype[so$celltype %in% c(3)] <- 'NK cell'
so$celltype[so$celltype %in% c(0, 4, 7)] <- 'Macrophage'
so$celltype[so$celltype %in% c(2, 5)] <- 'Monocyte'
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
<img width="840" height="840" alt="output_11_4" src="https://github.com/user-attachments/assets/604ae2c7-39da-4916-ab6e-dffbefe49e70" />
<img width="840" height="840" alt="output_11_3" src="https://github.com/user-attachments/assets/dbcaaa74-9603-4a24-9d12-c5f95d0185ae" />
<img width="840" height="840" alt="output_11_1" src="https://github.com/user-attachments/assets/a1dc38a9-167e-4b66-8550-91a0c0da98af" />
<img width="840" height="840" alt="output_11_0" src="https://github.com/user-attachments/assets/7c1ad12e-b1c3-4c46-9515-443374974284" />


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
<img width="840" height="840" alt="output_14_1" src="https://github.com/user-attachments/assets/22a4d5de-249b-4af5-b507-f7f0e2e63527" />
<img width="840" height="840" alt="output_14_0" src="https://github.com/user-attachments/assets/838275d4-8cae-4490-9849-ba2894d7a3cb" />

