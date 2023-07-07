from typing import List, Tuple, Optional
import numpy as np
import scanpy as sc
import squidpy as sq
import pandas as pd

def filter_flow_vars(adata: sc.AnnData,
                    vars_subset: List[str],
                    flowsig_expr_key: str = 'X_flow',
                    flowsig_network_key: str = 'flowsig_network'):
    
    X_flow_orig = adata.obsm[flowsig_expr_key]
    flowsig_info_orig = adata.uns[flowsig_network_key]

    # We define the list in this weird way to preserve all the var types etc
    flow_vars_to_subset = [flow_var for flow_var in flowsig_info['flow_vars'] if flow_var in vars_subset]
    subset_indices = [flowsig_info['flow_vars'].index(flow_var) for flow_var in flowsig_info['flow_vars'] if flow_var in vars_subset]

    X_flow = X_flow_orig[:, subset_indices]

    # Subset the flowsig network info as well
    flowsig_info = {}

    flow_variable_types = {}
    flow_downstream_tfs = {}
    flow_interactions = {}

    for flow_var in flow_vars_to_subset:

        flow_variable_types[flow_var] = flowsig_info_orig['flow_var_types'][flow_var]
        flow_downstream_tfs[flow_var] = flowsig_info_orig['downstream_tfs'][flow_var]
        flow_interactions[flow_var] = flowsig_info_orig['interactions'][flow_var]

    # Store the new FlowSig variable information
    flowsig_info = {'flow_vars': flow_vars_to_subset,
                        'flow_var_types': flow_variable_types,
                        'downstream_tfs': flow_downstream_tfs,
                        'interactions': flow_interactions
                        }

    adata.obsm[flowsig_expr_key + '_orig'] = X_flow_orig
    adata.obsm[flowsig_expr_key] = X_flow
    adata.uns[flowsig_network_key + '_orig'] = flowsig_info_orig
    adata.uns[flowsig_network_key] = flowsig_info

def determine_differentially_flowing_vars(adata: sc.AnnData,
                                        condition_key: str,
                                        control_key: str,
                                        flowsig_expr_key: str = 'X_flow',
                                        flowsig_network_key: str = 'flowsig_network',
                                        qval_threshold: float = 0.05,
                                        logfc_threshold: float = 0.5):
    
    # Construct AnnData for flow expression
    conditions = adata.obs[condition_key].unique().tolist()
    adata_flow = sc.AnnData(X=adata.obsm[flowsig_expr_key])
    adata_flow.var_names = pd.Index(adata.uns[flowsig_network_key]['flow_vars'])
    adata_flow.obs = adata.obs

    adata_flow.uns['log1p'] = {'base':None} # Just in case
    sc.tl.rank_genes_groups(adata_flow, key_added=condition_key, groupby=condition_key, method='wilcoxon')

    # Determine the differentially flowing vars
    differentially_flowing_vars = []

    lowqval_des = {}
    for cond in conditions:

        group_key = control_key + '_' + cond

        # Get the DEs with respect to this contrast
        result = sc.get.rank_genes_groups_df(adata_flow, group=group_key, key=condition_key).copy()
        result["-logQ"] = -np.log(result["pvals"].astype("float"))
        lowqval_de = result.loc[(abs(result["logfoldchanges"]) > logfc_threshold)&(abs(result["pvals_adj"]) < qval_threshold)]

        lowqval_des[cond] = lowqval_de['names'].tolist()
        
    differentially_flowing_vars = list(set.union(*map(set, [lowqval_des[cond] for cond in lowqval_des])))

    filter_flow_vars(adata,
                    differentially_flowing_vars,
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
                                    moran_threshold: float = 0.1,
                                    qval_threshold: float = 0.05,
                                    logfc_threshold: float = 0.5,
                                    coord_type: str = 'grid',
                                    n_neighbours: int = 6,
                                    library_key: str = None):
    

    if spatial: # We calculate the spatial autocorrelation (using Moran's I) and cut off genes below a defined threshold

        determine_spatially_flowing_vars(adata,
                                         flowsig_expr_key,
                                         flowsig_network_key,
                                         moran_threshold,
                                         coord_type,
                                         n_neighbours,
                                         library_key)
        
    else:

        determine_differentially_flowing_vars(adata, flowsig_expr_key,
                                            condition_key,
                                            qval_threshold, logfc_threshold)