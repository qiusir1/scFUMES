import scanpy as sc
import pandas as pd


def prepare_adata(
    adata,
    cluster_label="cell_type",
    min_genes=200,
    min_cells=10,
    log1p=True,
    regress_covariates: list[str] | None = None,
):
    """
    Prepare AnnData object for scFUMES analysis.
    """
    #if adata is already preprocessed
    if adata.raw is None:
        adata.raw = adata.copy()

    # QC & log‑normalisation
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.normalize_total(adata, target_sum=1e4)
    if log1p:   
        sc.pp.log1p(adata)
    if regress_covariates:
        missing_covariates = [
            covariate for covariate in regress_covariates if covariate not in adata.obs.columns
        ]
        if missing_covariates:
            missing = ", ".join(missing_covariates)
            raise ValueError(
                "The following covariates were not found in adata.obs: "
                f"{missing}. Ensure the columns exist before regression."
            )
        sc.pp.regress_out(adata, regress_covariates)
    adata.var_names_make_unique()
    
    print(f"Data contains {adata.n_obs} cells and {adata.n_vars} genes after preprocessing.")
    print(adata)
    
    return adata
