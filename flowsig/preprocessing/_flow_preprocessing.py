from __future__ import annotations
from typing import Tuple, Optional, Sequence
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

    X = adata.obsm[flowsig_network_key]
    if var_mask is not None:
        X = X[:, var_mask]

    adata_flow = AnnData(X=X, obs=adata.obs.copy())
    adata_flow.var = adata.uns[flowsig_network_key]["flow_var_info"].copy()
    if var_mask is not None:
        adata_flow.var = adata_flow.var.iloc[var_mask]
    return adata_flow

def _select_vars_union(
    rank_df_per_group: dict[str, pd.DataFrame],
    logfc_thr: float,
    qval_thr: float,
) -> list[str]:
    """Return union of gene names passing both thresholds in any group."""
    selected: set[str] = set()
    for df in rank_df_per_group.values():
        keep = (df["pvals_adj"] < qval_thr) & (np.abs(df["logfoldchanges"]) > logfc_thr)
        selected |= set(df.loc[keep, "names"]) # this is equivalent to selected = selected | set(df.loc[keep, "names"])
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
    logfc_thr: float = 0.5,
    qval_thr: float = 0.05,
    method: str = 'v1'
) -> None:
    
    flowsig_network_key = config.flowsig_expr_key
    if flowsig_network_key not in adata.uns:
        raise ValueError(f"Could not find {flowsig_network_key} in adata.uns")

    pert_conds = [cond for cond in adata.obs[condition_key].unique() if cond != control]
    info = adata.uns[config.flowsig_network_key]["flow_var_info"]

    inflow_mask = info["Type"].eq("inflow").to_numpy()
    outflow_mask = info["Type"].eq("outflow").to_numpy()

    ad_in  = _make_flow_adata(adata, config, inflow_mask)
    ad_out = _make_flow_adata(adata, config, outflow_mask)

    for ad in (ad_in, ad_out):
        ad.uns["log1p"] = {"base": None}  # prevent scanpy warning
        sc.tl.rank_genes_groups(ad, groupby=condition_key, method="wilcoxon")

    def _collect(ad) -> list[str]:
        pert_cond = {cond: sc.get.rank_genes_groups_df(ad, group=cond) for cond in pert_conds}
        return _select_vars_union(pert_cond, logfc_thr, qval_thr)

    vars_keep = _collect(ad_in) + _collect(ad_out) + info[info["Type"] == "module"].index.tolist()
    filter_flow_vars(adata, vars_keep, config)

def determine_spatially_flowing_vars(
    adata: AnnData,
    *,
    config: FlowSigConfig = FlowSigConfig(),
    moran_thr: float = 0.1,
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
        spatially_varying_vars.extend(ad.uns["moranI"].query("I > @moran_thr").index)

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