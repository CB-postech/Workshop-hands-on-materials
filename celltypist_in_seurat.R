library(sceasy)
library(reticulate)

# please use normalized Seurat object
celltypist_in_seurat <- function(so, conda_env, save_path, model_path) {
    h5ad_path = file.path(save_path, paste0("/tmp_adata.h5ad"))

    use_condaenv(conda_env)
    loompy <- reticulate::import('loompy')

    so[["RNA"]] <- as(so[["RNA"]], "Assay") # to avoid sceasy error occured by format conversion

    sceasy::convertFormat(so, from="seurat", to="anndata", drop_single_values=FALSE, outFile= h5ad_path)

    py_run_string("import scanpy as sc")
    py_run_string("import pandas as pd")
    py_run_string("import celltypist")
    py_run_string("from celltypist import models")
    py_run_string("from scipy import sparse")

    py_run_string(paste0(
        "adata = sc.read_h5ad('",
        h5ad_path,
        "')"
    ))

    py_run_string("sc.pp.normalize_per_cell(
                adata, counts_per_cell_after=10**4)")
    py_run_string("sc.pp.log1p(adata)")
    py_run_string(paste0("adata.X = sparse.csr_matrix(adata.X)"))

    # run celltypist
    py_run_string(
    paste0(
        "celltypist_result = celltypist.annotate(",
        "adata, majority_voting=True, over_clustering='seurat_clusters', ",
        "model='", model_path, "'",
        ")"
    )
    )

    # plot celltypist dotplots and save results
    py_code <- paste0(
        "celltypist.dotplot(celltypist_result,use_as_reference='seurat_clusters',use_as_prediction='majority_voting',save='celltypist_dotplot_majority_voting.png')\n",
        "celltypist.dotplot(celltypist_result,use_as_reference='seurat_clusters',use_as_prediction='predicted_labels',save='celltypist_dotplot_predicted_labels.png')\n",
        "adata_predicted = celltypist_result.to_adata()\n",
        "adata_predicted.obs[['predicted_labels','over_clustering','majority_voting','conf_score']].to_csv('",save_path,"/celltypist_result.csv',index=True)\n"
    )
    py_run_string(py_code)
}
