import scanpy as sc
import pandas as pd
import flowsig as fs

adata = sc.read('data/chen22_svg_E9.5.h5ad')

# We construct 10 gene expression modules using the raw cell count.
fs.construct_gems_using_nsf(adata,
                            n_gems = 20,
                            layer_key = 'count',
                            length_scale = 5.0)

commot_output_key = 'commot-cellchat'

# For spatial data, we need to construct spatial blocks that are used for block bootstrapping, 
# to preserve the spatial correlation of the gene expression data. The idea is that by sampling
# within these spatial blocks, we will better preserve these spatial correlation structures
# during bootstrapping. We construct the blocks using simple K-Means clustering over the spatial 
# locations.
fs.construct_spatial_blocks(adata,
                             n_blocks=20,
                             use_graph=False,
                             spatial_block_key = "spatial_block",
                             spatial_key = "spatial")

# We first construct the potential cellular flows from the commot output
fs.construct_flows_from_commot(adata,
                                commot_output_key,
                                gem_expr_key = 'X_gem',
                                scale_gem_expr = True,
                                flowsig_network_key = 'flowsig_network',
                                flowsig_expr_key = 'X_flow')

# Then we subset for "spatially flowing" variables
# (How do we turn the squidpy arguments into a kwargs)
fs.determine_informative_variables(adata,  
                                    flowsig_expr_key = 'X_flow',
                                    flowsig_network_key = 'flowsig_network',
                                    spatial = True,
                                    moran_threshold = 0.15,
                                    coord_type = 'grid',
                                    n_neighbours = 8,
                                    library_key = None)

# For spatial data, we need to construct spatial blocks that are used
# for block bootstrapping, to preserve the spatial correlation of the
# gene expression data.


# Now we are ready to learn the network
fs.learn_intercellular_flows(adata,
                        flowsig_key = 'flowsig_network',
                        flow_expr_key = 'X_flow',
                        use_spatial = True,
                        block_key = 'spatial_block',
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

adata.write('data/chen22_svg_E9.5.h5ad', compression='gzip')