from typing import List, Tuple, Optional
import numpy as np
import scanpy as sc
import squidpy as sq
import pandas as pd

def subset_for_flow_type(adata: sc.AnnData,
                         var_type: str = 'all',
                         flowsig_expr_key: str = 'X_flow',
                         flowsig_network_key: str = 'flowsig_network'):
    
    var_types = ['all', 'inflow', 'module', 'outflow']

    if var_type not in var_types:
        ValueError("Need to specify var_type as one of the following: %s"  % var_types)

    X_flow = adata.obsm[flowsig_expr_key]
    adata_subset = sc.AnnData(X=X_flow)
    adata_subset.obs = adata.obs
    adata_subset.var = adata.uns[flowsig_network_key]['flow_var_info']

    if var_type != 'all':

        adata_subset = adata_subset[:, adata_subset.var['Type'] == var_type]

    return adata_subset


def filter_flow_vars(adata: sc.AnnData,
                    vars_subset: List[str],
                    flowsig_expr_key: str = 'X_flow',
                    flowsig_network_key: str = 'flowsig_network'):
    
    X_flow_orig = adata.obsm[flowsig_expr_key]
    flow_var_info_orig = adata.uns[flowsig_network_key]['flow_var_info']
    flowsig_vars_orig = flow_var_info_orig.index.tolist()
    
    # We define the list in this weird way to preserve all the var types etc
    subset_indices = [flowsig_vars_orig.index(flow_var) for flow_var in vars_subset]

    X_flow = X_flow_orig[:, subset_indices]

    # Subset the flowsig network info as well
    flow_var_info = flow_var_info_orig[flow_var_info_orig.index.isin(vars_subset)]

    # Store the new FlowSig variable information
    flowsig_info = {'flow_var_info': flow_var_info}

    adata.obsm[flowsig_expr_key + '_orig'] = X_flow_orig
    adata.obsm[flowsig_expr_key] = X_flow
    adata.uns[flowsig_network_key + '_orig'] = flow_var_info_orig
    adata.uns[flowsig_network_key] = flowsig_info


def determine_differentially_flowing_vars(adata: sc.AnnData,
                                        condition_key: str,
                                        control_key: str,
                                        flowsig_expr_key: str = 'X_flow',
                                        flowsig_network_key: str = 'flowsig_network',
                                        qval_threshold: float = 0.05,
                                        logfc_threshold: float = 0.5):
    
    # Construct AnnData for flow expression
    perturbed_conditions = [cond for cond in adata.obs[condition_key].unique().tolist() if cond != control_key]

    # Construct inflow and outflow adata objects
    adata_inflow = subset_for_flow_type(adata,
                                        var_type = 'inflow',
                                        flowsig_expr_key = flowsig_expr_key,
                                        flowsig_network_key = flowsig_network_key)

    adata_outflow = subset_for_flow_type(adata,
                                        var_type = 'outflow',
                                        flowsig_expr_key = flowsig_expr_key,
                                        flowsig_network_key = flowsig_network_key)
       
    # Calculate differentially inflowing vars
    adata_inflow.uns['log1p'] = {'base':None} # Just in case
    sc.tl.rank_genes_groups(adata_inflow, key_added=condition_key, groupby=condition_key, method='wilcoxon')

    # Determine the differentially flowing vars
    diff_inflow_vars = []

    lowqval_des_inflow = {}
    for cond in perturbed_conditions:

        # Get the DEs with respect to this contrast
        result = sc.get.rank_genes_groups_df(adata_inflow, group=cond, key=condition_key).copy()
        result["-logQ"] = -np.log(result["pvals"].astype("float"))
        lowqval_de = result.loc[(abs(result["logfoldchanges"]) > logfc_threshold)&(abs(result["pvals_adj"]) < qval_threshold)]

        lowqval_des_inflow[cond] = lowqval_de['names'].tolist()
        
    diff_inflow_vars = list(set.union(*map(set, [lowqval_des_inflow[cond] for cond in lowqval_des])))

    # Calculate differentially inflowing vars
    adata_outflow.uns['log1p'] = {'base':None} # Just in case
    sc.tl.rank_genes_groups(adata_outflow, key_added=condition_key, groupby=condition_key, method='wilcoxon')

    # Determine the differentially flowing vars
    diff_outflow_vars = []

    lowqval_des_outflow = {}
    for cond in perturbed_conditions:

        # Get the DEs with respect to this contrast
        result = sc.get.rank_genes_groups_df(adata_outflow, group=cond, key=condition_key).copy()
        result["-logQ"] = -np.log(result["pvals"].astype("float"))
        lowqval_de = result.loc[(abs(result["logfoldchanges"]) > logfc_threshold)&(abs(result["pvals_adj"]) < qval_threshold)]

        lowqval_des_outflow[cond] = lowqval_de['names'].tolist()
        
    diff_outflow_vars = list(set.union(*map(set, [lowqval_des_outflow[cond] for cond in lowqval_des_outflow])))

    # We don't change GEM vars because from experience, they typically incorporate condition-specific changes as is
    gem_vars = 
    vars_to_subset = diff_inflow_vars + diff_outflow_vars + gem_vars

    filter_flow_vars(adata,
                    vars_to_subset,
                    flowsig_expr_key,
                    flowsig_network_key)
    
def determine_spatially_flowing_vars(adata: sc.AnnData,
                                    flowsig_expr_key: str = 'X_flow',
                                    flowsig_network_key: str = 'flowsig_network',
                                    moran_threshold: float = 0.1,
                                    coord_type: str = 'grid',
                                    n_neighbours: int = 6,
                                    library_key: str = None):

    
    # Construct AnnData for flow expression
    adata_flow = sc.AnnData(X=adata.obsm[flowsig_expr_key])
    adata_flow.var_names = pd.Index(adata.uns[flowsig_network_key]['flow_vars'])
    
    if 'spatial' not in adata.obsm:
        ValueError("Need to specify spatial coordinates in adata.obsm['spatial'].")
    else:
        adata_flow.obsm['spatial'] = adata.obsm['spatial']

        # Can't have spatial connectivities without spatial coordinates, lol
        if 'spatial_connectivities' not in adata.obsp['spatial_connectivities']:

            coord_types = ['grid', 'generic']
            
            if coord_type not in coord_types:
                ValueError("Please specify coord_type to be one of %s" % coord_types)

            sq.spatial_neighbors(adata_flow, coord_type=coord_type, n_neighs=n_neighbours, library_key=library_key)

            # Filter genes based on moran_threshold
            spatially_flowing_vars = adata_flow.uns['moranI_ligand'][adata_flow.uns['moranI_ligand']['I'] > moran_threshold].index.tolist()

            # Re-adjust the flow variables
            filter_flow_vars(adata,
                             spatially_flowing_vars,
                             flowsig_expr_key,
                             flowsig_network_key)

def determine_informative_variables(adata: sc.AnnData,  
                                    flowsig_expr_key: str = 'X_flow',
                                    flowsig_network_key: str = 'flowsig_network',
                                    spatial: bool = False,
                                    condition_key: str = None,
                                    control_key: str = None,
                                    moran_threshold: float = 0.1,
                                    qval_threshold: float = 0.05,
                                    logfc_threshold: float = 0.5,
                                    coord_type: str = 'grid',
                                    n_neighbours: int = 6,
                                    library_key: str = None):
    

    if spatial: # We calculate the spatial autocorrelation (using Moran's I) and cut off genes below a defined threshold

        determine_spatially_flowing_vars(adata,
                                         flowsig_expr_key=flowsig_expr_key,
                                         flowsig_network_key=flowsig_network_key,
                                         moran_threshold=moran_threshold,
                                         coord_type=coord_type,
                                         n_neighbours=n_neighbours,
                                         library_key=library_key)
        
    else:

        determine_differentially_flowing_vars(adata,
                                              condition_key=condition_key,
                                              control_key=control_key,
                                              flowsig_expr_key=flowsig_expr_key,
                                              flowsig_network_key=flowsig_network_key,
                                              qval_threshold=qval_threshold,
                                              logfc_threshold=logfc_threshold)
