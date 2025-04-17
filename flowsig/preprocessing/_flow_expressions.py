from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Literal, Any, Tuple
import numpy as np
import scanpy as sc
from anndata import AnnData
import pandas as pd
import os
from functools import lru_cache
import pathlib
import logging
logger = logging.getLogger(__name__)

@dataclass(frozen=True)
class FlowSigConfig:
    gem_expr_key: str = 'X_gem'
    scale_gem_expr: bool = True
    flowsig_network_key: str = 'flowsig_network'
    flowsig_expr_key: str = 'X_flow'

@lru_cache(maxsize=2)
def _load_cellchat_tfs(model_organism: str) -> pd.DataFrame:
    path = pathlib.Path(__file__).with_name(
        f'../data/cellchat_interactions_tfs_{model_organism}.csv.gz'
    )
    return pd.read_csv(path, index_col=0)

def _safe_get(adata: AnnData, gene: str) -> Optional[np.ndarray]:
    try:
        return adata[:, gene].X.A1          
    except KeyError:                        
        return None

def _dense_expr(adata: AnnData, genes: List[str]|str) -> np.ndarray:

    if isinstance(genes, str):
        return adata[:, genes].X.A1
    else:
        return adata[:, genes].X.toarray()      
    
def _assemble_flows(
        adata: AnnData,
        adata_outflow: AnnData,
        adata_inflow: AnnData,
        adata_gem: AnnData,
        *,
        config: FlowSigConfig
) -> None:
    # horizontal concat (dense or sparse)
    flow_expressions = np.hstack([
        adata_outflow.X,
        adata_inflow.X,
        adata_gem.X
    ])

    flow_variables = (
        adata_outflow.var_names.tolist()
        + adata_inflow.var_names.tolist()
        + adata_gem.var_names.tolist()
    )

    flow_var_info = pd.DataFrame({
        'Type': adata_outflow.var['type']
                  .append(adata_inflow.var['type'])
                  .append(adata_gem.var['type']),
        'Downstream_TF': adata_outflow.var['downstream_tfs']
                           .append(adata_inflow.var['downstream_tfs'])
                           .append(adata_gem.var['downstream_tfs']),
        'Interaction': adata_outflow.var['interactions']
                         .append(adata_inflow.var['interactions'])
                         .append(adata_gem.var['interactions']),
    }, index=pd.Index(flow_variables))

    adata.obsm[config.flowsig_expr_key] = flow_expressions
    adata.uns[config.flowsig_network_key] = {'flow_var_info': flow_var_info}

def construct_gem_expressions(adata: AnnData,
                            gem_expr_key: str = 'X_gem',
                            scale_gem_expr: bool = True,
                            layer_key: Optional[str] = None) -> Tuple[AnnData, List[str]]:
    
    gem_expressions = adata.obsm[gem_expr_key]

    # Scale so that the GEM memberships sum to 1 per cell
    gem_sum = gem_expressions.sum(axis=0)
    gem_expressions = gem_expressions / gem_sum
    
    num_gems = gem_expressions.shape[1]
    flow_gems = ['GEM-' + str(i + 1) for i in range(num_gems)]

    adata_gem = AnnData(X=gem_expressions)
    adata_gem.var.index = pd.Index(flow_gems)
    adata_gem.var['downstream_tfs'] = '' # For housekeeping for later
    adata_gem.var['type'] = 'module' # Define variable types
    adata_gem.var['interactions'] = '' # For housekeeping for later

    if scale_gem_expr:

        if layer_key is not None:

            scale_factor = adata.layers[layer_key].copy().sum(1).mean()

        else:
            scale_factor = np.expm1(adata.X).sum(1).mean()

        adata_gem.X *= scale_factor
        sc.pp.log1p(adata_gem)

    return adata_gem, flow_gems

def construct_inflow_signals_cellchat(adata: AnnData,
                                    cellchat_output_key: str, 
                                    model_organism: Literal['human', 'mouse'] = 'human',
                                    tfs_to_use: Optional[List[str]] = None,
                                    method: Literal['v1', 'v2'] = 'v1') -> Tuple[AnnData, List[str]]:

    vars_set = set(adata.var_names)  

    cellchat_interactions_and_tfs = _load_cellchat_tfs(model_organism)

    if cellchat_output_key not in adata.uns:
        raise KeyError(f"'{cellchat_output_key}' not found in adata.uns")

    ccc_output_merged = pd.concat([adata.uns[cellchat_output_key][sample] for sample in adata.uns[cellchat_output_key]])
    ccc_interactions = ccc_output_merged['interaction_name_2'].unique().tolist()

    unique_inflow_vars_and_interactions = {}

    # The inflow variables are constructed from receptor gene expression
    # that will be weighted by the average expression of their downstream
    # TF targets. These can be multi-units, but because of the weighting,
    # we hestitate to call them receptors in the context of cellular flows.
    for i, interaction in enumerate(ccc_interactions):
        receptor = interaction.split(' - ')[1].strip('()')
        receptor_split = receptor.split('+')
        receptors = []
        
        for i, rec in enumerate(receptor_split):
            if rec not in vars_set:
                receptor_v2_split = cellchat_interactions_and_tfs[cellchat_interactions_and_tfs['interaction_name_2'] == interaction]['receptor.symbol'].unique()[0].split(', ')
                
                receptors.append(receptor_v2_split[i])
            else:
                receptors.append(rec)
        
        receptor = '+'.join(receptors)

        if receptor not in unique_inflow_vars_and_interactions:
            unique_inflow_vars_and_interactions[receptor] = [interaction]
            
        else:
            unique_inflow_vars_and_interactions[receptor].append(interaction)
            
    inflow_vars = sorted(list(unique_inflow_vars_and_interactions.keys()))
    split_receptors = [receptor.split('+') for receptor in inflow_vars]
    unique_receptors = sorted({unit for rec in split_receptors for unit in rec})

    receptor_expression = _dense_expr(adata, unique_receptors)

    receptor_indices = {rec: i for i, rec in enumerate(unique_receptors)}
    # Take the log to speed up when we take the geometric mean
    log_receptor_expr = np.log(receptor_expression + 1e-12)

    inflow_expressions = np.empty((adata.n_obs, len(inflow_vars)))
    for k, units in enumerate(split_receptors):
        
        unit_cols = [receptor_indices[unit] for unit in units]

        # More efficient way of takingg geometric mean        
        inflow_expressions[:, k] = np.exp(log_receptor_expr[:, unit_cols].mean(axis=1))
        
    inflow_expressions_adjusted = inflow_expressions.copy()

    inflow_interactions = []
    unique_inflow_vars_and_tfs = {}

    for i, receptor in enumerate(inflow_vars):
        
        interactions_for_receptor = unique_inflow_vars_and_interactions[receptor]
        
        relevant_interactions = cellchat_interactions_and_tfs[cellchat_interactions_and_tfs['interaction_name_2'].isin(interactions_for_receptor)]
        
        joined_interactions = '/'.join(sorted(interactions_for_receptor))
        inflow_interactions.append(joined_interactions)
        
        downstream_tfs = []
        
        # Go through each category
        possible_downstream_tfs = relevant_interactions['Receptor-TF-combined'].dropna().tolist()\
                                    +  relevant_interactions['Ligand-TF-combined'].dropna().tolist()
        
        for unit in possible_downstream_tfs:
            split_unit = unit.split('_')
            for tf in split_unit:
                if tf not in downstream_tfs:
                    downstream_tfs.append(tf)
                    
        downstream_tfs = np.intersect1d(downstream_tfs, adata.var_names)

        if tfs_to_use is not None:
            downstream_tfs = np.intersect1d(downstream_tfs, tfs_to_use)
        
        unique_inflow_vars_and_tfs[receptor] = sorted(list(downstream_tfs))
        
        if len(downstream_tfs) != 0:
                            
            average_tf_expression = _dense_expr(adata, downstream_tfs).mean(axis=1).ravel()
                            
            inflow_expressions_adjusted[:, i] *= average_tf_expression
            
    inflow_downstream_tfs = []
    for i, inflow_var in enumerate(inflow_vars):
        
        interactions_for_receptor = unique_inflow_vars_and_interactions[inflow_var]
        downstream_tfs = unique_inflow_vars_and_tfs[inflow_var]
        inflow_downstream_tfs.append('_'.join(sorted(downstream_tfs)))
        
    adata_inflow = AnnData(X=inflow_expressions_adjusted)
    adata_inflow.var.index = pd.Index(inflow_vars)
    adata_inflow.var['downstream_tfs'] = inflow_downstream_tfs
    adata_inflow.var['type'] = 'inflow' # Define variable types
    adata_inflow.var['interactions'] = inflow_interactions

    return adata_inflow, inflow_vars

def construct_outflow_signals_cellchat(adata: AnnData,
                                    cellchat_output_key: str, 
                                    ) -> Tuple[AnnData, List[str]]:
    
    cellchat_output_merged = pd.concat([adata.uns[cellchat_output_key][sample] for sample in adata.uns[cellchat_output_key]])
    cellchat_interactions = cellchat_output_merged['interaction_name_2'].unique().tolist()
    relevant_interactions = {}

    for inter in cellchat_interactions:

        ligand = inter.split(' - ')[0]

        if ligand not in relevant_interactions:

            ligand_expr = _safe_get(adata, ligand)

            # Check if the alternative ligand name (sometimes CellChat is inconsistent)
            if ligand_expr is None:

                ligand = cellchat_output_merged[cellchat_output_merged['interaction_name_2'] == inter]['ligand'].values[0]

                ligand_expr = _safe_get(adata, ligand)
                    
            if ligand_expr is not None and ligand not in relevant_interactions:
                relevant_interactions[ligand] = [inter]

        else:
            relevant_interactions[ligand].append(inter)

    outflow_vars = list(relevant_interactions.keys())

    outflow_expressions = _dense_expr(adata, outflow_vars)

    adata_outflow = AnnData(X=outflow_expressions)
    adata_outflow.var.index = pd.Index(outflow_vars)
    adata_outflow.var['downstream_tfs'] = '' # For housekeeping for later
    adata_outflow.var['type'] = 'outflow' # Define variable types

    relevant_interactions_of_ligands = []
    for ligand in outflow_vars:
        interactions_of_ligand = relevant_interactions[ligand]
        relevant_interactions_of_ligands.append('/'.join(interactions_of_ligand))
        
    adata_outflow.var['interactions'] = relevant_interactions_of_ligands # Define variable types

    return adata_outflow, outflow_vars

def construct_flows_from_cellchat(adata: AnnData,
                                cellchat_output_key: str,
                                model_organism: Literal['human', 'mouse'] = 'human',
                                tfs_to_use: Optional[List[str]] = None,
                                config: FlowSigConfig = FlowSigConfig(),
                                method: Literal['v1', 'v2'] = 'v1') -> None:

    model_organisms = ['human', 'mouse']

    if model_organism not in model_organisms:
        raise ValueError ("Invalid model organism. Please select one of: %s" % model_organisms)
    
    if cellchat_output_key not in adata.uns:
        raise KeyError(f"'{cellchat_output_key}' not found in adata.uns")
    
    # Define the expression
    adata_outflow, outflow_vars = construct_outflow_signals_cellchat(adata, cellchat_output_key)

    adata_inflow, inflow_vars = construct_inflow_signals_cellchat(adata, cellchat_output_key, model_organism, tfs_to_use, method)

    adata_gem, flow_gem_vars = construct_gem_expressions(adata, config.gem_expr_key, config.scale_gem_expr)

    _assemble_flows(adata, adata_outflow, adata_inflow, adata_gem, config=config)
    
def construct_inflow_signals_cellphonedb(adata: AnnData,
                                    cellphonedb_output_key: str,
                                    cellphonedb_active_tfs_key: str) -> Tuple[AnnData, List[str]]:
    
    vars_set = set(adata.var_names)

    ccc_output_merged = pd.concat([adata.uns[cellphonedb_output_key][sample] for sample in adata.uns[cellphonedb_output_key]])
    cpdb_active_tfs_merged = pd.concat([adata.uns[cellphonedb_active_tfs_key][sample] for sample in adata.uns[cellphonedb_active_tfs_key]])
    ccc_interactions = ccc_output_merged['interacting_pair'].unique().tolist()

    unique_inflow_vars_and_interactions = {}

    # The inflow variables are constructed from receptor gene expression
    # that will be weighted by the average expression of their downstream
    # TF targets. These can be multi-units, but because of the weighting,
    # we hestitate to call them receptors in the context of cellular flows.
    for i, interaction in enumerate(ccc_interactions):
        receptors = ccc_output_merged[ccc_output_merged['interacting_pair'] == interaction]['gene_b'].dropna().unique().tolist()
        
        for receptor in receptors:
            if receptor in vars_set:
                if receptor not in unique_inflow_vars_and_interactions:
                    unique_inflow_vars_and_interactions[receptor] = [interaction]
                    
                else:
                    unique_inflow_vars_and_interactions[receptor].append(interaction)
                
    inflow_vars = sorted(list(unique_inflow_vars_and_interactions.keys()))
            
    inflow_expressions = _dense_expr(adata, inflow_vars)

    inflow_expressions_adjusted = inflow_expressions.copy()

    inflow_interactions = []
    unique_inflow_vars_and_tfs = {}

    for i, receptor in enumerate(inflow_vars):
        
        interactions_for_receptor = unique_inflow_vars_and_interactions[receptor]
        inflow_interactions.append('/'.join(sorted(interactions_for_receptor)))
        
        # Go through each category
        possible_downstream_tfs = cpdb_active_tfs_merged[cpdb_active_tfs_merged['gene_b'] == receptor]['active_TF'].dropna().unique().tolist()
                    
        if len(possible_downstream_tfs) != 0:
            downstream_tfs = np.intersect1d(possible_downstream_tfs, adata.var_names)
        
            unique_inflow_vars_and_tfs[receptor] = sorted(list(downstream_tfs))

        else:
            unique_inflow_vars_and_tfs[receptor] = ''
        
        if len(possible_downstream_tfs) != 0:
            
            average_tf_expression = _dense_expr(adata, possible_downstream_tfs).mean(1).ravel()
                            
            inflow_expressions_adjusted[:, i] *= average_tf_expression
            
    inflow_downstream_tfs = []
    for i, inflow_var in enumerate(inflow_vars):
        
        interactions_for_receptor = unique_inflow_vars_and_interactions[inflow_var]
        downstream_tfs = unique_inflow_vars_and_tfs[inflow_var]
        inflow_downstream_tfs.append('_'.join(sorted(downstream_tfs)))
        
    adata_inflow = AnnData(X=inflow_expressions_adjusted)
    adata_inflow.var.index = pd.Index(inflow_vars)
    adata_inflow.var['downstream_tfs'] = inflow_downstream_tfs
    adata_inflow.var['type'] = 'inflow' # Define variable types
    adata_inflow.var['interactions'] = inflow_interactions

    return adata_inflow, inflow_vars

def construct_outflow_signals_cellphonedb(adata: AnnData,
                                    cellphonedb_output_key: str) -> Tuple[AnnData, List[str]]:

    vars_set = set(adata.var_names)

    cellphonedb_output_merged = pd.concat([adata.uns[cellphonedb_output_key][sample] for sample in adata.uns[cellphonedb_output_key]])
    cellphonedb_interactions = cellphonedb_output_merged['interacting_pair'].unique().tolist()
    relevant_interactions = {}

    for inter in cellphonedb_interactions:

        ligands = cellphonedb_output_merged[cellphonedb_output_merged['interacting_pair'] == inter]['gene_a'].dropna().unique().tolist()

        for ligand in ligands:
            if (ligand in vars_set):

                if (ligand not in relevant_interactions):
                    relevant_interactions[ligand] = [inter]

                else:
                    relevant_interactions[ligand].append(inter)

    outflow_vars = list(relevant_interactions.keys())

    outflow_expressions = _dense_expr(adata, outflow_vars)

    adata_outflow = AnnData(X=outflow_expressions)
    adata_outflow.var.index = pd.Index(outflow_vars)
    adata_outflow.var['downstream_tfs'] = '' # For housekeeping for later
    adata_outflow.var['type'] = 'outflow' # Define variable types

    relevant_interactions_of_ligands = []
    for ligand in outflow_vars:
        interactions_of_ligand = relevant_interactions[ligand]
        relevant_interactions_of_ligands.append('/'.join(interactions_of_ligand))
        
    adata_outflow.var['interactions'] = relevant_interactions_of_ligands # Define variable types

    return adata_outflow, outflow_vars

def construct_flows_from_cellphonedb(adata: AnnData,
                                cellphonedb_output_key: str,
                                cellphonedb_tfs_key: str,
                                model_organism: str = 'human',
                                config: FlowSigConfig = FlowSigConfig()):

    if model_organism != 'human':
        ValueError("CellPhoneDB only supports human data.")

    if cellphonedb_output_key not in adata.uns:
        raise KeyError(f"'{cellphonedb_output_key}' not found in adata.uns")
    
    if cellphonedb_tfs_key not in adata.uns:
        raise KeyError(f"'{cellphonedb_tfs_key}' not found in adata.uns")   

    # Define the expression
    adata_outflow, outflow_vars = construct_outflow_signals_cellphonedb(adata, cellphonedb_output_key)

    adata_inflow, inflow_vars = construct_inflow_signals_cellphonedb(adata, cellphonedb_output_key, cellphonedb_tfs_key)

    adata_gem, flow_gem_vars = construct_gem_expressions(adata, config.gem_expr_key, config.scale_gem_expr)

    _assemble_flows(adata, adata_outflow, adata_inflow, adata_gem, config=config)
    
def construct_inflow_signals_liana(adata: AnnData,
                                    liana_output_key: str, 
                                    use_tfs: bool = False,
                                    model_organism: Literal['human', 'mouse'] = 'human') -> Tuple[AnnData, List[str]]:
    
    model_organisms = ['human', 'mouse']

    if model_organism not in model_organisms:
        raise ValueError ("Invalid model organism. Please select one of: %s" % model_organisms)

    vars_set = set(adata.var_names)
    
    if use_tfs:

        cellchat_interactions_and_tfs = _load_cellchat_tfs(model_organism)

        ccc_output_merged = pd.concat([adata.uns[liana_output_key][sample] for sample in adata.uns[liana_output_key]])
        inflow_vars = sorted(ccc_output_merged['receptor_complex'].unique().tolist())
        inflow_vars = [var.replace('_', '+') for var in inflow_vars]

        unique_inflow_vars_and_interactions = {}

        # The inflow variables are constructed from receptor gene expression
        # that will be weighted by the average expression of their downstream
        # TF targets. These can be multi-units, but because of the weighting,
        # we hestitate to call them receptors in the context of cellular flows.
        for inflow_var in inflow_vars:
            relevant_rec = inflow_var.replace('+', '_')

            relevant_ligands = sorted(ccc_output_merged[ccc_output_merged['receptor_complex'] == relevant_rec]['ligand_complex'].unique().tolist())
            interactions_for_receptor = []

            for lig in relevant_ligands:
                if '+' in inflow_var:
                    interactions_for_receptor.append(lig + ' - ' + '(' + inflow_var + ')')
                else:
                    interactions_for_receptor.append(lig + ' - ' + inflow_var)

            unique_inflow_vars_and_interactions[inflow_var] = interactions_for_receptor
                
        inflow_vars = sorted(list(unique_inflow_vars_and_interactions.keys()))
        split_receptors = [receptor.split('+') for receptor in inflow_vars]
        unique_receptors = sorted({unit for rec in split_receptors for unit in rec})

        receptor_expression = _dense_expr(adata, unique_receptors)

        receptor_indices = {rec: i for i, rec in enumerate(unique_receptors)}
        # Take the log to speed up when we take the geometric mean
        log_receptor_expr = np.log(receptor_expression + 1e-12)

        inflow_expressions = np.empty((adata.n_obs, len(inflow_vars)))
        for k, units in enumerate(split_receptors):
            
            unit_cols = [receptor_indices[unit] for unit in units]

            # More efficient way of takingg geometric mean        
            inflow_expressions[:, k] = np.exp(log_receptor_expr[:, unit_cols].mean(axis=1))
        
        inflow_expressions_adjusted = inflow_expressions.copy()

        inflow_interactions = []
        unique_inflow_vars_and_tfs = {}

        for i, receptor in enumerate(inflow_vars):
            
            interactions_for_receptor = unique_inflow_vars_and_interactions[receptor]
            
            relevant_interactions = cellchat_interactions_and_tfs[cellchat_interactions_and_tfs['interaction_name_2'].isin(interactions_for_receptor)]
            
            joined_interactions = '/'.join(sorted(interactions_for_receptor))
            inflow_interactions.append(joined_interactions)
            
            downstream_tfs = []
            
            # Go through each category
            possible_downstream_tfs = relevant_interactions['Receptor-TF-combined'].dropna().tolist()
            
            for unit in possible_downstream_tfs:
                split_unit = unit.split('_')
                for tf in split_unit:
                    if tf not in downstream_tfs:
                        downstream_tfs.append(tf)
                        
            downstream_tfs = np.intersect1d(downstream_tfs, adata.var_names)
            
            unique_inflow_vars_and_tfs[receptor] = sorted(list(downstream_tfs))
            
            if len(downstream_tfs) != 0:
                
                average_tf_expression = _dense_expr(adata, downstream_tfs).mean(axis=1).ravel()
                                    
                inflow_expressions_adjusted[:, i] *= average_tf_expression
                
        inflow_downstream_tfs = []
        for i, inflow_var in enumerate(inflow_vars):
            
            interactions_for_receptor = unique_inflow_vars_and_interactions[inflow_var]
            downstream_tfs = unique_inflow_vars_and_tfs[inflow_var]
            inflow_downstream_tfs.append('_'.join(sorted(downstream_tfs)))
            
        adata_inflow = AnnData(X=inflow_expressions_adjusted)
        adata_inflow.var.index = pd.Index(inflow_vars)
        adata_inflow.var['downstream_tfs'] = inflow_downstream_tfs
        adata_inflow.var['type'] = 'inflow' # Define variable types
        adata_inflow.var['interactions'] = inflow_interactions

        return adata_inflow, inflow_vars
    
def construct_outflow_signals_liana(adata: AnnData,
                                    liana_output_key: str, 
                                    ) -> Tuple[AnnData, List[str]]:
    vars_set = set(adata.var_names)

    liana_output_merged = pd.concat([adata.uns[liana_output_key][sample] for sample in adata.uns[liana_output_key]])
   
    outflow_vars = sorted(liana_output_merged['ligand_complex'].unique().tolist())
    outflow_vars = [var.replace('_', '+') for var in outflow_vars]

    relevant_interactions = {}

    for outflow_var in outflow_vars:
            
        relevant_lig = outflow_var.replace('+', '_')

        interactions_with_ligand = []

        relevant_receptors = sorted(liana_output_merged[liana_output_merged['ligand_complex'] == relevant_lig]['receptor_complex'].unique().tolist())

        for rec in relevant_receptors:
            
            inter = outflow_var + ' - ' + rec.replace('_', '+')
            interactions_with_ligand.append(inter)

        relevant_interactions[outflow_var] = interactions_with_ligand

    outflow_expressions = np.zeros((adata.n_obs, len(outflow_vars)))
    split_ligands = [ligand.split('+') for receptor in outflow_vars]
    unique_ligands = sorted({unit for lig in split_ligands for unit in lig})

    ligand_expression = _dense_expr(adata, unique_ligands)

    ligand_indices = {lig: i for i, lig in enumerate(unique_ligands)}
    # Take the log to speed up when we take the geometric mean
    log_ligand_expr = np.log(ligand_expression + 1e-12)

    outflow_expressions = np.empty((adata.n_obs, len(outflow_vars)))
    for k, units in enumerate(split_ligands):
        
        unit_cols = [ligand_indices[unit] for unit in units]

        # More efficient way of takingg geometric mean        
        outflow_expressions[:, k] = np.exp(log_ligand_expr[:, unit_cols].mean(axis=1))

    adata_outflow = AnnData(X=outflow_expressions)
    adata_outflow.var.index = pd.Index(outflow_vars)
    adata_outflow.var['downstream_tfs'] = '' # For housekeeping for later
    adata_outflow.var['type'] = 'outflow' # Define variable types

    relevant_interactions_of_ligands = []
    for ligand in outflow_vars:
        interactions_of_ligand = relevant_interactions[ligand]
        relevant_interactions_of_ligands.append('/'.join(interactions_of_ligand))
        
    adata_outflow.var['interactions'] = relevant_interactions_of_ligands # Define variable types

    return adata_outflow, outflow_vars
   
def construct_flows_from_liana(adata: AnnData,
                                liana_output_key: str,
                                use_tfs: bool = False,
                                model_organism: Literal['human', 'mouse'] = 'human',
                                config: FlowSigConfig = FlowSigConfig()) -> None:

    if liana_output_key not in adata.uns:
        raise KeyError(f"'{liana_output_key}' not found in adata.uns")
    
    # Define the expression
    adata_outflow, outflow_vars = construct_outflow_signals_liana(adata, liana_output_key)

    adata_inflow, inflow_vars = construct_inflow_signals_liana(adata, liana_output_key, use_tfs, model_organism)

    adata_gem, flow_gem_vars = construct_gem_expressions(adata, config.gem_expr_key, config.scale_gem_expr)

    _assemble_flows(adata, adata_outflow, adata_inflow, adata_gem, config=config)

def construct_inflow_signals_commot(adata: AnnData,
                                    commot_output_key: str) -> Tuple[AnnData, List[str]]:

    # Inflow variables are inferred from outflow variables
    outflow_vars = sorted(adata.uns[commot_output_key + '-info']['df_ligrec']['ligand'].unique().tolist())
    inflow_vars = ['inflow-' + outflow_var for outflow_var in outflow_vars]


    inflow_interactions = []
    inflow_expressions = np.empty((adata.n_obs, len(inflow_vars)))
    for i, inflow_var in enumerate(inflow_vars):
        lig = inflow_var.strip('inflow-')
        inferred_interactions = [pair.replace('r-', '') for pair in adata.obsm[commot_output_key + '-sum-receiver'].columns if pair.startswith('r-' + lig)]
        inflow_interactions.append('/'.join(sorted(inferred_interactions)))

        # We sum the total received signal across each interaction at each spot
        cols = [f"r-{inter}" for inter in inferred_interactions]
        inflow_expressions[:, i] = adata.obsm[f"{commot_output_key}-sum-receiver"][cols].sum(1)

    adata_inflow = AnnData(X=inflow_expressions)
    adata_inflow.var.index = pd.Index(inflow_vars)
    adata_inflow.var['downstream_tfs'] = ''
    adata_inflow.var['type'] = 'inflow' 
    adata_inflow.var['interactions'] = inflow_interactions

    return adata_inflow, inflow_vars

def construct_outflow_signals_commot(adata: AnnData,
                                    commot_output_key: str
                                    ) -> Tuple[AnnData, List[str]]:
    
    vars_set = set(adata.var_names)
    # Inflow variables are inferred from outflow variables
    outflow_vars = sorted(adata.uns[commot_output_key + '-info']['df_ligrec']['ligand'].unique().tolist())
    outflow_vars = [var for var in outflow_vars if var in vars_set]

    outflow_interactions = []

    for i, outflow_var in enumerate(outflow_vars):

        inferred_interactions = [pair[2:] for pair in adata.obsm[commot_output_key + '-sum-sender'].columns if pair.startswith('s-' + outflow_var)]
        outflow_interactions.append('/'.join(sorted(inferred_interactions)))

    # Outflow signal expression is simply ligand gene expression
    outflow_expressions = _dense_expr(adata, outflow_vars)

    adata_outflow = AnnData(X=outflow_expressions)
    adata_outflow.var.index = pd.Index(outflow_vars)
    adata_outflow.var['downstream_tfs'] = ''
    adata_outflow.var['type'] = 'outflow' 
    adata_outflow.var['interactions'] = outflow_interactions

    return adata_outflow, outflow_vars

def construct_flows_from_commot(adata: AnnData,
                                commot_output_key: str,
                                config: FlowSigConfig = FlowSigConfig()) -> None:

    if commot_output_key not in adata.uns:
        raise KeyError(f"'{commot_output_key}' not found in adata.uns")
    
    # Define the expression
    adata_outflow, outflow_vars = construct_outflow_signals_commot(adata, commot_output_key)

    adata_inflow, inflow_vars = construct_inflow_signals_commot(adata, commot_output_key)

    adata_gem, flow_gem_vars = construct_gem_expressions(adata, config.gem_expr_key, config.scale_gem_expr)

    # Determine the flow_variables
    _assemble_flows(adata, adata_outflow, adata_inflow, adata_gem, config=config)

def construct_flow_expressions(
        adata: AnnData,
        *,
        spatial: bool,
        method: Literal['cellchat', 'cellphonedb', 'liana'] | None = None,
        **kwargs: Any,
) -> None:

    if spatial:
        # For spatial data, flows are based on COMMOT output
        construct_flows_from_commot(adata, **kwargs)
        return

    else:
        if method is None:
            raise ValueError(
                "Argument 'method' must be provided when spatial=False. "
                "Choose from: 'cellchat', 'cellphonedb', or 'liana'."
            )
        
        method = method.lower()
        if method == 'cellchat':
            construct_flows_from_cellchat(adata, **kwargs)
        elif method == 'cellphonedb':
            construct_flows_from_cellphonedb(adata, **kwargs)
        elif method == 'liana':
            construct_flows_from_liana(adata, **kwargs)
        else:
            raise ValueError(
                f"Unrecognised method '{method}'. "
                "Valid options are: 'cellchat', 'cellphonedb', or 'liana'."
            )