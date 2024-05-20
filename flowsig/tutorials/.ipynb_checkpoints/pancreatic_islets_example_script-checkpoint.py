import scanpy as sc
import pandas as pd
import flowsig as fs

adata = sc.read('burkhardt21_merged.h5ad')
condition_key = 'Condition'

# We construct 10 gene expression modules using the raw cell count.
fs.construct_gems_using_pyliger(adata,
                                n_gems = 10,
                                layer_key = 'counts',
                                condition_key = condition_key)

# Make sure your keys for these align with their condition labels
cellchat_Ctrl = pd.read('burkhardt21_leiden_communications_Ctrl.csv')
cellchat_IFNg = pd.read('burkhardt21_leiden_communications_IFNg.csv')

cellchat_output_key = 'cellchat_output'
adata.uns[cellchat_output_key] = {'Ctrl': cellchat_Ctrl,
                                  'IFNg': cellchat_IFNg}

# We first construct the potential cellular flows from the cellchat output
fs.construct_flows_from_cellchat(adata,
                                cellchat_output_key,
                                gem_expr_key = 'X_gem',
                                scale_gem_expr = True,
                                model_organism = 'human',
                                flowsig_network_key = 'flowsig_network',
                                flowsig_expr_key = 'X_flow')

# Then we subset for "differentially flowing" variables
fs.determine_informative_variables(adata,  
                                    flowsig_expr_key = 'X_flow',
                                    flowsig_network_key = 'flowsig_network',
                                    spatial = False,
                                    condition_key = condition_key,
                                    qval_threshold = 0.05,
                                    logfc_threshold = 0.5)

# Now we are ready to learn the network
fs.learn_intercellular_flows(adata,
                        condition_key = condition_key,
                        control_key = 'Ctrl', 
                        flowsig_key = 'flowsig_network',
                        flow_expr_key = 'X_flow',
                        use_spatial = False,
                        n_jobs = 1,
                        n_bootstraps = 10)

# Now we do post-learning validation to reorient the network and remove low-quality edges.
# This part is key for reducing false positives
fs.apply_biological_flow(adata,
                        flowsig_network_key = 'flowsig_network',
                        adjacency_key = 'adjacency',
                        validated_adjacency_key = 'adjacency_validated')

edge_threshold = 0.7

fs.filter_low_confidence_edges(adata,
                                edge_threshold = edge_threshold,
                                flowsig_network_key = 'flowsig_network',
                                adjacency_key = 'adjacency',
                                filtered_adjacency_key = 'adjacency_filtered')

adata.write('data/burkhardt21_merged.h5ad', compression='gzip')