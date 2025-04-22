from __future__ import annotations
from typing import Tuple, Optional, Sequence, Literal
import numpy as np
import scanpy as sc
from anndata import AnnData
import squidpy as sq
import pandas as pd
from ..preprocessing import FlowSigConfig

def _make_flow_adata(
    adata: AnnData,
    config: FlowSigConfig,
    var_mask: np.ndarray | None = None,
) -> AnnData:

    flowsig_expr_key = config.flowsig_expr_key
    flowsig_network_key = config.flowsig_network_key

    if flowsig_expr_key not in adata.obsm:
        raise ValueError(f"Could not find {flowsig_expr_key} in adata.obsm")
    
    if flowsig_network_key not in adata.uns:
        raise ValueError(f"Could not find {flowsig_network_key} in adata.uns")

    X = adata.obsm[flowsig_expr_key]

    adata_flow = AnnData(X=X, obs=adata.obs.copy())
    adata_flow.var = adata.uns[flowsig_network_key]["flow_var_info"].copy()

    if var_mask is not None:
        adata_flow = adata_flow[:, var_mask]
        
    return adata_flow

def _log_fold_cohens(adata: AnnData,
                     condition_key: str,
                     group1: str,
                     group2: str) -> pd.DataFrame:
    """Calculate a standardised log-fold change somewhat equivalent to Cohen's d between group1 and group2.
    This is adapted from the scoreMarkers function in scran (see here: https://rdrr.io/github/MarioniLab/scran/man/scoreMarkers.html).
    """

    # Calculate the mean and standard deviation for each group
    expr_group1 = adata[adata.obs[condition_key] == group1].X
    expr_group2 = adata[adata.obs[condition_key] == group2].X
    mean1 = expr_group1.mean(axis=0)
    mean2 = expr_group2.mean(axis=0)
    std1 = expr_group1.std(axis=0)
    std2 = expr_group2.std(axis=0)

    log_means_diff = mean1 - mean2
    lfc_cohen = log_means_diff / ( 0.5 * (std1**2.0 + std2**2.0) )** 0.5
    lfc_cohens_results = pd.DataFrame(data=lfc_cohen, index=adata.var_names, columns=['logfoldchanges_cohen'])

    return lfc_cohens_results

def _select_vars_union(
    logfoldchanges_per_group: dict[str, pd.DataFrame],
    logfc_threshold: float,
    qval_threshold: float = None,
    construction: Literal['v1', 'v2'] = 'v1'
) -> list[str]:
    """Return union of gene names passing both thresholds in any group."""
    selected: set[str] = set()

    if construction == 'v1':
        for df in logfoldchanges_per_group.values():
            keep = (np.abs(df["logfoldchanges"]) > logfc_threshold) & (df["pvals_adj"] < qval_threshold)

            selected |= set(df.loc[keep, "names"]) # this is equivalent to selected = selected | set(df.loc[keep, "names"])

    else:
        for df in logfoldchanges_per_group.values():
            keep = (np.abs(df["logfoldchanges_cohen"]) > logfc_threshold)

            selected |= set(df.index[keep])

    return list(selected)

def subset_for_flow_type(
    adata: AnnData,
    var_type: str = "all",
    config: FlowSigConfig = FlowSigConfig(),
) -> AnnData:

    flowsig_network_key = config.flowsig_network_key

    if flowsig_network_key not in adata.uns:
        raise ValueError(f"Could not find {flowsig_network_key} in adata.uns")
    
    if var_type not in {"all", "inflow", "module", "outflow"}:
        raise ValueError("var_type must be 'all', 'inflow', 'module', or 'outflow'.")

    info = adata.uns[config.flowsig_network_key]["flow_var_info"]
    mask = None if var_type == "all" else (info["Type"].to_numpy() == var_type)
    return _make_flow_adata(adata, config, mask)

def filter_flow_vars(
    adata: AnnData,
    vars_subset: Sequence[str],
    config: FlowSigConfig = FlowSigConfig(),
) -> None:
    
    flowsig_network_key = config.flowsig_network_key
    flowsig_expr_key = config.flowsig_expr_key

    if flowsig_network_key not in adata.uns:
        raise ValueError(f"Could not find {flowsig_network_key} in adata.uns")
    if flowsig_expr_key not in adata.obsm:  
        raise ValueError(f"Could not find {flowsig_expr_key} in adata.obsm")
    
    full_info = adata.uns[config.flowsig_network_key]["flow_var_info"]
    mask = full_info.index.isin(vars_subset)

    # no change → return early
    if mask.all():
        return

    # backup originals once
    adata.obsm.setdefault(f"{config.flowsig_expr_key}_orig", adata.obsm[config.flowsig_expr_key])
    adata.uns.setdefault(f"{config.flowsig_network_key}_orig", {"flow_var_info": full_info})

    adata.obsm[flowsig_expr_key] = adata.obsm[config.flowsig_expr_key][:, mask]
    adata.uns[flowsig_network_key] = {"flow_var_info": full_info.loc[mask]}

def determine_differentially_flowing_vars(
    adata: AnnData,
    condition_key: str,
    control: str,
    *,
    config: FlowSigConfig = FlowSigConfig(),
    logfc_threshold: float = None,
    qval_threshold: float = None,
    construction: Literal['v1', 'v2'] = 'v1'
) -> None:
    
    flowsig_expr_key = config.flowsig_expr_key
    if flowsig_expr_key not in adata.obsm:
        raise ValueError(f"Could not find {flowsig_expr_key} in adata.obsm")

    flowsig_network_key = config.flowsig_network_key
    if flowsig_network_key not in adata.uns:
        raise ValueError(f"Could not find {flowsig_network_key} in adata.uns")

    pert_conds = [cond for cond in adata.obs[condition_key].unique() if cond != control]
    info = adata.uns[config.flowsig_network_key]["flow_var_info"]

    inflow_mask = info["Type"].eq("inflow").to_numpy()
    outflow_mask = info["Type"].eq("outflow").to_numpy()

    ad_in  = _make_flow_adata(adata, config, inflow_mask)
    ad_out = _make_flow_adata(adata, config, outflow_mask)

    vars_keep = None

    def _collect(ad, construction) -> list[str]:
            
        pert_cond = {}
        if construction == 'v1':
            pert_cond = {cond: sc.get.rank_genes_groups_df(ad, group=cond) for cond in pert_conds}
        
        else: # Just select based on logfc_thr for logfoldchanges_cohen
            pert_cond = {cond: _log_fold_cohens(ad, condition_key=condition_key, group1=cond, group2=control) for cond in pert_conds}

        return _select_vars_union(pert_cond, logfc_threshold, qval_threshold, construction)

    if construction == 'v1':

        for ad in (ad_in, ad_out):
            ad.uns["log1p"] = {"base": None}  # prevent scanpy warning
            sc.tl.rank_genes_groups(ad, groupby=condition_key, method="wilcoxon")

    vars_keep = _collect(ad_in, construction) + _collect(ad_out, construction) + info[info["Type"] == "module"].index.tolist()


    filter_flow_vars(adata, vars_keep, config)

def determine_spatially_flowing_vars(
    adata: AnnData,
    *,
    config: FlowSigConfig = FlowSigConfig(),
    moran_threshold: float = 0.1,
    coord_type: str = "grid",
    n_neighs: int = 6,
    library_key: str | None = None,
    n_perms: int | None = None,
    n_jobs: int | None = None,
) -> None:
    
    if "spatial" not in adata.obsm:
        raise ValueError("adata.obsm['spatial'] not found.")

    flowsig_network_key = config.flowsig_network_key
    if flowsig_network_key not in adata.uns:
        raise ValueError(f"Could not find {flowsig_network_key} in adata.uns")
    
    info = adata.uns[flowsig_network_key]["flow_var_info"]
    masks = {var_type: info["Type"].eq(var_type).to_numpy() for var_type in ("inflow", "outflow")}

    adatas = {
        var_type: _make_flow_adata(adata, config, mask) for var_type, mask in masks.items()
    }
    for ad in adatas.values():
        ad.obsm["spatial"] = adata.obsm["spatial"]

    # build neighbours only once per coord_type/param set
    for ad in adatas.values():
        if "spatial_connectivities" not in ad.obsp:
            sq.gr.spatial_neighbors(ad, coord_type=coord_type, n_neighs=n_neighs, library_key=library_key)
        sq.gr.spatial_autocorr(ad, genes=ad.var_names, n_perms=n_perms, n_jobs=n_jobs)

    spatially_varying_vars = []
    for var_type, ad in adatas.items():
        spatially_varying_vars.extend(ad.uns["moranI"].query("I > @moran_threshold").index)

    vars_keep = spatially_varying_vars + info[info["Type"] == "module"].index.tolist()
    filter_flow_vars(adata, vars_keep, config)

def determine_informative_variables(
    adata: AnnData,
    *,
    config: FlowSigConfig = FlowSigConfig(),
    spatial: bool = False,
    condition_key: str | None = None,
    control: str | None = None,
    **kwargs,
) -> None:
    
    """Wrapper that chooses spatial vs. differential selection."""
    if spatial:
        determine_spatially_flowing_vars(adata, config=config, **kwargs)
    else:
        if condition_key is None or control is None:
            raise ValueError("condition_key and control must be provided for differential selection.")
        determine_differentially_flowing_vars(
            adata,
            condition_key,
            control,
            config=config,
            **kwargs,
        )