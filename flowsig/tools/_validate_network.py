from typing import Optional, Dict
import scanpy as sc
import networkx as nx
import numpy as np
import pandas as pd
import graphical_models as gpm
import networkx as nx 

def filter_low_confidence_edges(adata: sc.AnnData,
                                edge_threshold: float,
                                flowsig_network_key: str = 'flowsig_network',
                                adjacency_key: str = 'adjacency',
                                filtered_key: str = 'filtered'):
    """
    Validate the learned CPDAG from UT-IGSP by checking edges against the assumed
    biological flow model, inflow signal -> gene expression module -> outflow signal.
    As the CPDAG contains directed arcs and undirected edges, we remove directed arcs
    that 

    Parameters
    ----------
    adata
        The annotated dataframe (typically from Scanpy) of the single-cell data.
        You need to have run FlowSig before running this step.

    edge_threshold
        The relative frequency of bootstrap edge frequency above which we keep edges.
        For directed arcs, we consider single edge frequencies. For undirected edges,
        we consider total edge weight.

    flowsig_network_key 
        The label in adata.uns where all of the flowsig output is stored, including the learned
        adjacency corresponding to the CPDAG (markov equivalence class), the flow variables used
        for inference, as well as their "flow types", which can be either inflow (ing), module
        (TFs or factors), or outflow (ing) signals.

    adjacency_key
        String key that specifies the adjacency for the learned network is stored in adata.uns[flowsig_network_key].

    adjacency_filtered_key
        String key that specifies where the validated network will be stored.

    Returns
    -------
    adjacency_filtered
        Matrix that encodes the CPDAG containing high-confidence directed arcs and 
        undirected arcs.

    """

    # Get the adjacency
    adjacency = adata.uns[flowsig_network_key]['network'][adjacency_key]
    flow_vars = adata.uns[flowsig_network_key]['flow_var_info'].index.tolist()
    flow_var_info = adata.uns[flowsig_network_key]['flow_var_info']

    cpdag = gpm.PDAG.from_amat(adjacency)

    adjacency_filtered = adjacency.copy()
    
    # First, let us calculate the total edge weights
    total_edge_weights = {}

    nonzero_rows, nonzero_cols = adjacency.nonzero()

    for i in range(len(nonzero_rows)):

        row_ind = nonzero_rows[i]
        col_ind = nonzero_cols[i]

        node_1 = flow_vars[row_ind]
        node_2 = flow_vars[col_ind]

        edge = (node_1, node_2)            

        if (edge not in total_edge_weights):
            total_edge_weights[edge] = adjacency[row_ind, col_ind]
        else:
            if (edge[1], edge[0]) in total_edge_weights:
                total_edge_weights[(edge[1], edge[0])] += adjacency[row_ind, col_ind]

    for arc in cpdag.arcs:

        node_1 = flow_vars[tuple(arc)[0]]
        node_2 = flow_vars[tuple(arc)[1]]

        edge_weight =  adjacency[row_ind, col_ind]

        # Need to account for both (node1, node2) and (node1, node1) as 
        # adjacency encodes directed network
        if edge_weight < edge_threshold: 
            adjacency_filtered[row_ind, col_ind] = 0.0

    for edge in cpdag.edges:

        node_1 = flow_vars[tuple(edge)[0]]
        node_2 = flow_vars[tuple(edge)[1]]

        # For directed arcs, we simply consider the total edge weights
        total_edge_weight = 0.0

        if (node_1, node_2) in total_edge_weights:

            total_edge_weight = total_edge_weights[(node_1, node_2)]

        else:

            total_edge_weight = total_edge_weights[(node_2, node_1)]

        # Need to account for both (node1, node2) and (node1, node1) as 
        # adjacency encodes directed network
        if total_edge_weight < edge_threshold: 

            adjacency_filtered[row_ind, col_ind] = 0.0
            adjacency_filtered[col_ind, row_ind] = 0.0
            
    # Save the "validated" adjacency
    filtered_adjacency_key = adjacency_key + '_' + filtered_key
    adata.uns[flowsig_network_key]['network'][filtered_adjacency_key] = adjacency_filtered

def apply_biological_flow(adata: sc.AnnData,
                        flowsig_network_key: str = 'flowsig_network',
                        adjacency_key: str = 'adjacency',
                        validated_key: str = 'validated'):
    """
    Validate the learned CPDAG from UT-IGSP by checking edges against the assumed
    biological flow model, inflow signal -> gene expression module -> outflow signal.
    As the CPDAG contains directed arcs and undirected edges, we remove directed arcs
    that 

    Parameters
    ----------
    adata
        The annotated dataframe (typically from Scanpy) of the single-cell data.
        You need to have run FlowSig before running this step.

    flowsig_network_key 
        The label in adata.uns where all of the flowsig output is stored, including the learned
        adjacency corresponding to the CPDAG (markov equivalence class), the flow variables used
        for inference, as well as their "flow types", which can be either inflow (ing), module (TFs or factors),
        or outflow (ing) signals.

    adjacency_key
        String key that specifies the adjacency for the learned network is stored in adata.uns[flowsig_network_key].

    adjacency_validated_key
        String key that specifies where the validated network will be stored.

    Returns
    -------
    adjacency_validated
        Matrix that encodes the CPDAG containing "biologically realistic" inferred flows, from
        inflow variables, to gene expression module variables, to outflow variables.

    """

    # Get the adjacency
    adjacency = adata.uns[flowsig_network_key]['network'][adjacency_key]
    flow_vars = adata.uns[flowsig_network_key]['flow_var_info'].index.tolist()
    flow_var_info = adata.uns[flowsig_network_key]['flow_var_info']

    cpdag = gpm.PDAG.from_amat(adjacency)
    adjacency_validated = adjacency.copy()
    
    for arc in cpdag.arcs:

        node_1 = flow_vars[tuple(arc)[0]]
        node_2 = flow_vars[tuple(arc)[1]]

        # Classify node types
        node_1_type = flow_var_info.loc[node_1]['Type']
        node_2_type = flow_var_info.loc[node_2]['Type']

        # Now we decide whether or not to add the  edges
        add_edge = False

        # Define the edge because we may need to reverse it
        edge = (node_1, node_2)

        # If there's a link from received morphogen to a TF
        if ( (node_1_type == 'inflow')&(node_2_type == 'module') ):

            add_edge = True

        if ( (node_1_type == 'module')&(node_2_type == 'outflow') ):

            add_edge = True

        if ((node_1_type == 'module')&(node_2_type == 'module')):

            add_edge = True

        if not add_edge:

            row_ind = list(flow_vars).index(node_1)
            col_ind = list(flow_vars).index(node_2)

            adjacency_validated[row_ind, col_ind] = 0.0

    for edge in cpdag.edges:

        node_1 = flow_vars[tuple(edge)[0]]
        node_2 = flow_vars[tuple(edge)[1]]

        # Classify node types
        node_1_type = flow_var_info.loc[node_1]['Type']
        node_2_type = flow_var_info.loc[node_2]['Type']

        # Define the edge because we may need to reverse it
        edge = (node_1, node_2)

        add_edge = False

        # If there's a link from received morphogen to a TF
        if ( (node_1_type == 'inflow')&(node_2_type == 'module') ):

            add_edge = True

        # If there's a link from received morphogen to a TF
        if ( (node_1_type == 'module')&(node_2_type == 'inflow') ):

            edge = (edge[1], edge[0])
            add_edge = True

        if ( (node_1_type == 'module')&(node_2_type == 'ligand') ):

            add_edge = True

        if ( (node_1_type == 'ligand')&(node_2_type == 'module') ):

            edge = (edge[1], edge[0])
            add_edge = True

        if ((node_1_type == 'module')&(node_2_type == 'module')):

            add_edge = True

        if not add_edge:

            row_ind = list(flow_vars).index(node_1)
            col_ind = list(flow_vars).index(node_2)

            # If we reversed the edge, make sure we get the right edge right:
            if edge[0] != node_1:

                row_ind = list(flow_vars).index(node_2)
                col_ind = list(flow_vars).index(node_1)

            adjacency_validated[row_ind, col_ind] = 0.0
    
    # Save the "validated" adjacency
    validated_adjacency_key = adjacency_key + '_' + validated_key
    adata.uns[flowsig_network_key]['network'][validated_adjacency_key] = adjacency_validated

def construct_intercellular_flow_network(adata: sc.AnnData,
                        flowsig_network_key: str = 'flowsig_network',
                        adjacency_key: str = 'adjacency_validated'):
    
    flow_vars = flow_vars
    flow_adjacency = adata.uns[flowsig_network_key]['network'][adjacency_key]

    nonzero_rows, nonzero_cols = flow_adjacency.nonzero()

    total_edge_weights = {}

    for i in range(len(nonzero_rows)):
        row_ind = nonzero_rows[i]
        col_ind = nonzero_cols[i]

        node_1 = flow_vars[row_ind]
        node_2 = flow_vars[col_ind]

        edge = (node_1, node_2)            

        if ( (edge not in total_edge_weights)&((edge[1], edge[0]) not in total_edge_weights) ):
            total_edge_weights[edge] = flow_adjacency[row_ind, col_ind]
        else:
            if (edge[1], edge[0]) in total_edge_weights:
                total_edge_weights[(edge[1], edge[0])] += flow_adjacency[row_ind, col_ind]
    
    
    flow_network = nx.DiGraph()

    # Now let's consturct the graph from the CPDAG
    cpdag =  gpm.PDAG.from_amat(flow_adjacency)

    # Add the directed edges (arcs) first
    for arc in cpdag.arcs:

        node_1 = flow_vars[tuple(arc)[0]]
        node_2 = flow_vars[tuple(arc)[1]]

        # Classify node types
        node_1_type = adata.uns[flowsig_network_key]['flow_var_types'][node_1]
        node_2_type = adata.uns[flowsig_network_key]['flow_var_types'][node_2]

        # Now we decide whether or not to add the damn edges
        add_edge = False

        # Define the edge because we may need to reverse it
        edge = (node_1, node_2)

        # If there's a link from the expressed morphogen to the received morphogen FOR the same morphogen
        if ( (node_1_type == 'inflow')&(node_2_type == 'module') ):
            add_edge = True

        # If there's a link from received morphogen to a TF
        if ( (node_1_type == 'module')&(node_2_type == 'outflow') ):

            add_edge = True

        if ( (node_1_type == 'module')&(node_2_type == 'module') ):

            add_edge = True

        if add_edge:

            # Get the total edge weight
            total_edge_weight = 0.0

            if edge in total_edge_weights:

                total_edge_weight = total_edge_weights[edge]

            else:

                total_edge_weight = total_edge_weights[(edge[1], edge[0])]

            edge_weight = flow_adjacency[tuple(arc)[0], tuple(arc)[1]]

            flow_network.add_edge(*edge)
            flow_network.edges[edge[0], edge[1]]['weight'] = edge_weight / total_edge_weight

    for edge in cpdag.edges:

        node_1 = flow_vars[tuple(arc)[0]]
        node_2 = flow_vars[tuple(arc)[1]]

        # Classify node types
        node_1_type = adata.uns[flowsig_network_key]['flow_var_types'][node_1]
        node_2_type = adata.uns[flowsig_network_key]['flow_var_types'][node_2]

        # Define the edge because we may need to reverse it
        edge = (node_1, node_2)

        add_edge = False
        # If there's a link from the expressed morphogen to the received morphogen FOR the same morphogen
        if ( (node_1_type == 'inflow')&(node_2_type == 'module') ):

            add_edge = True

        if ( (node_1_type == 'module')&(node_2_type == 'inflow') ):

            add_edge = True
            edge = (edge[1], edge[0])

        if ( (node_1_type == 'module')&(node_2_type == 'ligand') ):

            add_edge = True

        if ( (node_1_type == 'ligand')&(node_2_type == 'module') ):

            add_edge = True
            edge = (edge[1], edge[0])

        if ((node_1_type == 'module')&(node_2_type == 'module')):

            add_edge = True

        if add_edge:

            # Get the total edge weight
            total_edge_weight = 0.0

            if edge in total_edge_weights:

                total_edge_weight = total_edge_weights[edge]

            else:

                total_edge_weight = total_edge_weights[(edge[1], edge[0])]

            flow_network.add_edge(*edge)
            flow_network.edges[edge[0], edge[1]]['weight'] = min(total_edge_weight, 1.0)

            # Add the other way if we have modules
            if ((node_1_type == 'module')&(node_2_type == 'module')):

                flow_network.add_edge(edge[1], edge[0])
                flow_network.edges[edge[1], edge[0]]['weight'] = min(total_edge_weight, 1.0)

    return flow_network