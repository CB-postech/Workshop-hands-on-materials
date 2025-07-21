library(sceasy)
library(reticulate)

batch_aware_in_seurat <- function(so, batch_key = 'batch', nHVG = 2000, conda_env, save_path, is_normalized = TRUE) {
    h5ad_path = file.path(save_path, paste0("/tmp_adata.h5ad"))
    hvg_path = file.path(save_path, paste0("/hvg_", nHVG, "_", batch_key, ".csv"))

    use_condaenv(conda_env)
    loompy <- reticulate::import('loompy')

    so[["RNA"]] <- as(so[["RNA"]], "Assay") # to avoid sceasy error occured by format conversion

    sceasy::convertFormat(so, from="seurat", to="anndata", drop_single_values=FALSE, outFile= h5ad_path)

    py_run_string("import scanpy as sc")
    py_run_string("import pandas as pd")
    py_run_string(paste0(
        "adata = sc.read_h5ad('",
        h5ad_path,
        "')"
    ))
    if(!is_normalized) {
        py_run_string("sc.pp.normalize_total(adata, target_sum=1e4)")
        py_run_string("sc.pp.log1p(adata)")
    }

    py_run_string(paste0(
        "sc.pp.highly_variable_genes(adata, batch_key='",
        batch_key,
        "')"
    ))
    py_run_string(paste0(
        "pd.Series(adata.var_names[adata.var['highly_variable']], name='gene')",
        ".to_csv('",
        hvg_path,
        "', index=False)"
    ))
}
