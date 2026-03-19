import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import scanpy as sc
import numpy as np 
import pandas as pd
from typing import Union, Sequence, List, Tuple, Optional, Iterable
from matplotlib.axes import Axes

palette_network = list(sns.color_palette("tab20") \
                       + sns.color_palette("tab20b") \
                       + sns.color_palette("tab20c")\
                       + sns.color_palette("Set1")\
                       + sns.color_palette("Set2")\
                       + sns.color_palette('Set3')\
                       + sns.color_palette('Dark2')\
                       + sns.color_palette('Pastel1')\
                       + sns.color_palette('Pastel2'))

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
    adata_subset.var = pd.DataFrame(adata.uns[flowsig_network_key]['flow_var_info'])

    if var_type != 'all':

        adata_subset = adata_subset[:, adata_subset.var['Type'] == var_type]

    return adata_subset

def label_point(x, y, val, ax):

    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)

    for i, point in a.iterrows():

        ax.text(point['x']+0.1, point['y'] + 0.05, str(point['val']), fontdict={'size':12.0})


def plot_differentially_flowing_signals(adata: sc.AnnData,
                                        condition_key: str,
                                        pert_key: Optional[Iterable[str]] = None,
                                        var_type: str = 'all',
                                        flowsig_expr_key: str = 'X_flow',
                                        flowsig_network_key: str = 'flowsig_network',
                                        qval_threshold: float = 0.05,
                                        logfc_threshold: float = 0.5,
                                        label_lowqval: bool = False,
                                        scatter_size: float = 50,
                                        ax: Axes = None):

    var_types = ['all', 'inflow', 'module', 'outflow']

    adata_subset = subset_for_flow_type(adata = adata,
                                        var_type = var_type,
                                        flowsig_expr_key = flowsig_expr_key,
                                        flowsig_network_key = flowsig_network_key)
    
    if condition_key not in adata_subset.uns:
        adata_subset.uns['log1p'] = {'base': None}
        sc.tl.rank_genes_groups(adata_subset, key_added=condition_key, groupby=condition_key, method='wilcoxon')

    result = sc.get.rank_genes_groups_df(adata_subset, group=pert_key, key=condition_key).copy()
    # result["-logQ"] = -np.log(result["pvals_adj"].astype("float"))
    result["-logP"] = -np.log(result["pvals"].astype("float"))

    lowqval_de = result.loc[(abs(result["logfoldchanges"]) > logfc_threshold)&(abs(result["pvals_adj"]) < qval_threshold)]
    other_de = result.loc[(abs(result["logfoldchanges"]) <= logfc_threshold)|(abs(result["pvals_adj"]) >= qval_threshold)]

    if ax is None:
        fig, ax = plt.subplots()

    sns.regplot(
        x=other_de["logfoldchanges"],
        y=other_de["-logP"],
        fit_reg=False,
        scatter_kws={"s": scatter_size},
        ax=ax
    )
    sns.regplot(
        x=lowqval_de["logfoldchanges"],
        y=lowqval_de["-logP"],
        fit_reg=False,
        scatter_kws={"s": scatter_size},
        ax=ax
    )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)

    # Plot where the thresholds are
    ax.axhline(-np.log(qval_threshold), color='grey', alpha=0.5, linestyle='--')
    ax.axvline(-logfc_threshold, color='grey', alpha=0.5, linestyle='--')
    ax.axvline(logfc_threshold, color='grey', alpha=0.5, linestyle='--')

    ax.set_xlabel("log2 FC")
    ax.set_ylabel("-log p-value")

    if label_lowqval:
        label_point(lowqval_de['logfoldchanges'],
                    lowqval_de['-logP'],
                    lowqval_de['names'],
                    plt.gca())  
            
    ax.autoscale_view()
    return ax

def plot_intercellular_flows(adata: sc.AnnData,
                            flow_network: nx.DiGraph,
                            inflow_vars: Union[str, Sequence[str]] = None,
                            module_vars: Union[str, Sequence[str]] = None,
                            outflow_vars: Union[str, Sequence[str]] = None,
                            flowsig_network_key: str = 'flowsig_network',
                            align_mode: str = 'horizontal',
                            width_scale: float = 2.0,
                            x_margin_offset: float = 0.3,
                            y_margin_offset: float = 0.0,
                            ax: Axes = None):
    
    flow_vars = list(flow_network.nodes())
    flow_var_info = adata.uns[flowsig_network_key]['flow_var_info']

    inflows_used = inflow_vars is not None
    modules_used = module_vars is not None
    outflows_used = outflow_vars is not None

    var_types_used = np.array([inflows_used, modules_used, outflows_used], dtype=int)

    all_inflow_vars = [node for node in flow_vars if flow_var_info.loc[node]['Type'] == 'inflow']
    all_module_vars = [node for node in flow_vars if flow_var_info.loc[node]['Type'] == 'module']
    all_outflow_vars = [node for node in flow_vars if flow_var_info.loc[node]['Type'] == 'outflow']

    # The extremal cases are the easiest to consider
    if var_types_used.sum() == 0: # None of them have been considered so we plot all
        inflow_vars = all_inflow_vars
        module_vars = all_module_vars
        outflow_vars = all_outflow_vars

    # Only one was specified, so we have to trace the network to get the other two types
    elif var_types_used.sum() == 1:

        # If inflow vars have been specified, we get the GEMs connected to the specified inflows, and then the outflows connected to those inferred GEMs
        if inflows_used:
            module_vars = []
            for inflow in inflow_vars:
                inflow_successors = [node for node in flow_network.successors(inflow) if node in all_module_vars]

                for module in inflow_successors:
                    outflow_successors_of_module = [node for node in flow_network.successors(module) if node in all_outflow_vars]
                    if len(outflow_successors_of_module) != 0 and module not in module_vars:
                        module_vars.append(module)

            outflow_vars = []
            for module in module_vars:
                outflow_successors = [node for node in flow_network.successors(module) if node in all_outflow_vars]

                for outflow in outflow_successors:
                    if outflow not in outflow_vars:
                        outflow_vars.append(outflow)

        # If module vars have been specified, we get the inflow connected to the specified GEMs, and then the outflows connected to the specified GEMs
        elif modules_used:
            inflow_vars = []
            outflow_vars = []

            for module in module_vars:
                module_predecessors = [node for node in flow_network.predecessors(module) if node in all_inflow_vars]
                module_successors = [node for node in flow_network.successors(module) if node in all_outflow_vars]

                for inflow in module_predecessors:
                    if inflow not in inflow_vars:
                        inflow_vars.append(inflow)

                for outflow in module_successors:
                    if outflow not in outflow_vars:
                        outflow_vars.append(outflow)            
            
        # If outflow vars have been specified, we get the GEMs connected to the specified outflows, and then the inflow connected to the inferred GEMs
        else:
            module_vars = []
            for outflow in outflow_vars:
                outflow_predecessors = [node for node in flow_network.predecessors(outflow) if node in all_module_vars]

                for module in outflow_predecessors:
                    module_predecessors = [node for node in flow_network.predecessors(module) if node in all_inflow_vars]
                    if len(module_predecessors) != 0 and module not in module_vars:
                        module_vars.append(module)

            inflow_vars = []
            for module in module_vars:
                module_predecessors = [node for node in flow_network.predecessors(module) if node in all_inflow_vars]

                for inflow in module_predecessors:
                    if inflow not in inflow_vars:
                        inflow_vars.append(inflow)

    elif var_types_used.sum() == 2:

        # We just extract the outflows connected to the module vars
        if inflows_used and modules_used:
            outflow_vars = list(set([edge[1] for edge in flow_network.edges() if edge[0] in module_vars and edge[1] in all_outflow_vars]))
        
        # We need the union of GEMs connected to either the inflows or the outflows
        elif inflows_used and outflows_used:
            modules_from_inflows = list(set([edge[1] for edge in flow_network.edges() if edge[0] in inflow_vars and edge[1] in all_module_vars]))
            modules_to_outflows = list(set([edge[0] for edge in flow_network.edges() if edge[0] in all_module_vars and edge[1] in outflow_vars]))
            module_vars = list(set(modules_from_inflows + modules_to_outflows))
            
        # We get the inflows connected to the specified modules
        else: 
            inflow_vars = list(set([edge[0] for edge in flow_network.edges() if edge[0] in all_inflow_vars and edge[1] in module_vars]))

    # We will sort the inflow vars etc now
    compressed_inflows_used = all('-' in var for var in inflow_vars) # Assumption is that all inflow vars are written like 'Inflow-1', etc.
    compressed_modules_used = all('-' in var for var in module_vars) # Assumption is that all module vars are written like 'GEM-1', etc.
    compressed_outflows_used = all('-' in var for var in outflow_vars) # Assumption is that all outflow vars are written like 'Outflow-1', etc.
    inflow_vars = sorted(inflow_vars) if not compressed_inflows_used else sorted(inflow_vars, key=lambda x: int(x.split('-')[1])) # Assumption is that compressed inflow vars are written like 'Inflow-1', etc.
    module_vars = sorted(module_vars) if not compressed_modules_used else sorted(module_vars, key=lambda x: int(x.split('-')[1])) # Assumption is that all module vars are written like 'GEM-1', etc.
    outflow_vars = sorted(outflow_vars) if not compressed_outflows_used else sorted(outflow_vars, key=lambda x: int(x.split('-')[1])) # Assumption is that compressed outflow vars are written like 'Outflow-1', etc.

    resultant_pattern_graph = flow_network.subgraph(inflow_vars + module_vars + outflow_vars)

    node_colours = {}
    count = 0

    for inflow_var in all_inflow_vars:
        node_colours[inflow_var] = palette_network[count]
        count += 1

    for module_var in all_module_vars:
        node_colours[module_var] = palette_network[count]
        count += 1

    for ouflow_var in all_outflow_vars:
        node_colours[ouflow_var] = palette_network[count]
        count += 1

    for node in resultant_pattern_graph.nodes():

        if align_mode == 'horizontal':

            if node in inflow_vars:
                resultant_pattern_graph.nodes[node]['level'] = 3

            elif node in module_vars:
                resultant_pattern_graph.nodes[node]['level'] = 2

            else: # Downstream outflow vars
                resultant_pattern_graph.nodes[node]['level'] = 1

        else: # Otherwise we're going veritcally

            if node in inflow_vars:
                resultant_pattern_graph.nodes[node]['level'] = 1

            elif node in module_vars:
                resultant_pattern_graph.nodes[node]['level'] = 2

            else:
                resultant_pattern_graph.nodes[node]['level'] = 3

    # We will sort everything one last time out of paranoia
    resultant_pattern_graph_pos = nx.multipartite_layout(resultant_pattern_graph, subset_key='level', align=align_mode, scale=2.0)

    # Fix the damned positions
    inflow_pos = [resultant_pattern_graph_pos[node] for node in inflow_vars]
    module_pos = [resultant_pattern_graph_pos[node] for node in module_vars]
    outflow_pos = [resultant_pattern_graph_pos[node] for node in outflow_vars]

    if align_mode == 'horizontal':
        inflow_pos = sorted(inflow_pos, key=lambda x: x[0], reverse=False)
        module_pos = sorted(module_pos, key=lambda x: x[0], reverse=False)
        outflow_pos = sorted(outflow_pos, key=lambda x: x[0], reverse=False)

    else:
        inflow_pos = sorted(inflow_pos, key=lambda x: x[1], reverse=False)
        module_pos = sorted(module_pos, key=lambda x: x[1], reverse=False)
        outflow_pos = sorted(outflow_pos, key=lambda x: x[1], reverse=False)

    for i, node in enumerate(inflow_vars):
        resultant_pattern_graph_pos[node] = inflow_pos[i]
        
    for i, node in enumerate(module_vars):
        resultant_pattern_graph_pos[node] = module_pos[i]
        
    for i, node in enumerate(outflow_vars):
        resultant_pattern_graph_pos[node] = outflow_pos[i]

    # Let's sort the nodes by alphabetical order
    resultant_pattern_graph_edge_colours = [node_colours[edge[0]] for edge in resultant_pattern_graph.edges()]
    edge_widths = [width_scale*resultant_pattern_graph[u][v]['weight'] for u,v in resultant_pattern_graph.edges()]

    if ax is None:
        fig, ax = plt.subplots()
    nx.draw_networkx_edges(resultant_pattern_graph,
                            resultant_pattern_graph_pos,
                            edge_color=resultant_pattern_graph_edge_colours,
                            width=edge_widths,
                             alpha=0.75,
                             connectionstyle="arc3,rad=0.2",
                             ax=ax)
    nx.draw_networkx_nodes(resultant_pattern_graph,
                           resultant_pattern_graph_pos,
                           nodelist=inflow_vars, node_color=[node_colours[node] for node in inflow_vars], linewidths=1.0, edgecolors='black', node_size=100, ax=ax)
    nx.draw_networkx_nodes(resultant_pattern_graph,
                           resultant_pattern_graph_pos,
                           nodelist=module_vars, node_color=[node_colours[node] for node in module_vars], node_shape='v', node_size=100, ax=ax)
    nx.draw_networkx_nodes(resultant_pattern_graph,
                           resultant_pattern_graph_pos,
                           nodelist=outflow_vars, node_color=[node_colours[node] for node in outflow_vars], node_shape='s', node_size=100, ax=ax)

    inflow_labels = nx.draw_networkx_labels(resultant_pattern_graph,
                                            resultant_pattern_graph_pos,
                                            labels = {node:node for node in resultant_pattern_graph.nodes() if node in inflow_vars},
                                            font_size=14,
                                            font_family='Arial',
                                            ax=ax,
                                            horizontalalignment='left',
                                            verticalalignment = 'bottom');
    for text in inflow_labels.values():
        text.set_rotation(45)

    module_labels = nx.draw_networkx_labels(resultant_pattern_graph,
                                            resultant_pattern_graph_pos,
                                            labels = {node:node for node in resultant_pattern_graph.nodes() if node in module_vars},
                                            font_size=14, font_family='Arial',
                                            ax=ax,
                                            horizontalalignment='left',
                                            verticalalignment = 'bottom');
    for text in module_labels.values():
        text.set_rotation(45)

    outflow_labels = nx.draw_networkx_labels(resultant_pattern_graph,
                                            resultant_pattern_graph_pos,
                                            labels = {node:node for node in resultant_pattern_graph.nodes() if node in outflow_vars},
                                            font_size=14, font_family='Arial',
                                            ax=ax,
                                            horizontalalignment='right',
                                            verticalalignment = 'top');
    for text in outflow_labels.values():
        text.set_rotation(45)

    plt.tight_layout()
    plt.margins(x=x_margin_offset)
    plt.margins(y=y_margin_offset)
    plt.box(False)

    return ax

def plot_intercellular_cre_flows(adata: sc.AnnData,
                            flow_network: nx.DiGraph,
                            inflow_vars: Union[str, Sequence[str]] = None,
                            module_tf_vars: Union[str, Sequence[str]] = None,
                            module_cre_vars: Union[str, Sequence[str]] = None,
                            outflow_vars: Union[str, Sequence[str]] = None,
                            flowsig_network_key: str = 'flowsig_network',
                            align_mode: str = 'horizontal',
                            width_scale: float = 2.0,
                            x_margin_offset: float = 0.3,
                            y_margin_offset: float = 0.0,
                            ax: Axes = None):
    
    flow_vars = list(flow_network.nodes())
    flow_var_info = adata.uns[flowsig_network_key]['flow_var_info']

    inflows_used = inflow_vars is not None
    module_tfs_used = module_tf_vars is not None
    module_cres_used = module_cre_vars is not None
    outflows_used = outflow_vars is not None

    var_types_used = np.array([inflows_used, module_tfs_used, module_cres_used, outflows_used], dtype=int)

    all_inflow_vars = [node for node in flow_vars if flow_var_info.loc[node]['Type'] == 'inflow']
    all_module_tf_vars = [node for node in flow_vars if flow_var_info.loc[node]['Type'] == 'module_tf']
    all_module_cre_vars = [node for node in flow_vars if flow_var_info.loc[node]['Type'] == 'module_cre']
    all_outflow_vars = [node for node in flow_vars if flow_var_info.loc[node]['Type'] == 'outflow']

    # The extremal cases are the easiest to consider
    if var_types_used.sum() == 0: # None of them have been considered so we plot all
        inflow_vars = all_inflow_vars
        module_tf_vars = all_module_tf_vars
        module_cre_vars = all_module_cre_vars
        outflow_vars = all_outflow_vars

    # Only one was specified, so we have to trace the network to get the other three types
    elif var_types_used.sum() == 1:

        # If inflow vars have been specified, we get the TF modules connected to the specified inflows,
        # the CRE modules connected to those TF modules, and then the outflows connected to those inferred CRE modules
        if inflows_used:
            module_tf_vars = []
            module_cre_vars = []
            outflow_vars = []
            
            for inflow in inflow_vars:
                inflow_successors = [node for node in flow_network.successors(inflow) if node in all_module_tf_vars]

                for module_tf in inflow_successors:
                    module_cre_successors_of_module = [node for node in flow_network.successors(module_tf) if node in all_module_cre_vars]
                    if len(module_cre_successors_of_module) != 0 and module_tf not in module_tf_vars:
                        module_tf_vars.append(module_tf)

            for module_tf in module_tf_vars:
                module_cre_successors = [node for node in flow_network.successors(module_tf) if node in all_module_cre_vars]

                for module_cre in module_cre_successors:
                    outflow_successors_of_module_cre = [node for node in flow_network.successors(module_cre) if node in all_outflow_vars]
                    if len(outflow_successors_of_module_cre) != 0 and module_cre not in module_cre_vars:
                         module_cre_vars.append(module_cre) 

            for module_cre in module_cre_vars:
                outflow_successors = [node for node in flow_network.successors(module_cre) if node in all_outflow_vars]

                for outflow in outflow_successors:
                    if outflow not in outflow_vars:
                        outflow_vars.append(outflow)

        # If module vars have been specified, we get the inflow connected to the specified TF modules,
        # # the CRE modules connected to those TF modules, and then the outflows connected to the specified GEMs
        elif module_tfs_used:
            inflow_vars = []
            module_cre_vars = []
            outflow_vars = []

            for module_tf in module_tf_vars:
                module_tf_predecessors = [node for node in flow_network.predecessors(module_tf) if node in all_inflow_vars]
                module_cre_successors = [node for node in flow_network.successors(module_tf) if node in all_module_cre_vars]

                for inflow in module_tf_predecessors:
                    if inflow not in inflow_vars:
                        inflow_vars.append(inflow)

                for module_cre in module_cre_successors:
                    outflow_successor_of_module_cre = [node for node in flow_network.successors(module_cre) if node in all_outflow_vars]

                    if len(outflow_successors_of_module_cre) != 0 and module_cre not in module_cre_vars:
                        module_cre_vars.append(module_cre)
                            
                        for outflow in outflow_successor_of_module_cre:
                            if outflow not in outflow_vars:
                                outflow_vars.append(outflow)

        elif module_cres_used:
            inflow_vars = []
            module_tf_vars = []
            outflow_vars = []

            for module_cre in module_cre_vars:
                module_cre_predecessors = [node for node in flow_network.predecessors(module_cre) if node in all_module_tf_vars]
                module_cre_successors = [node for node in flow_network.successors(module_cre) if node in all_module_cre_vars]

                for module_tf in module_cre_predecessors:
                    inflow_predecessor_of_module_cre = [node for node in flow_network.predecessors(module_tf) if node in all_inflow_vars]

                    if len(inflow_predecessor_of_module_cre) != 0 and module_tf not in module_tf_vars:
                        module_tf_vars.append(module_tf)

                        for inflow in inflow_predecessor_of_module_cre:
                            if inflow not in inflow_vars:
                                inflow_vars.append(inflow)

                for module_cre in module_cre_successors:
                    outflow_successor_of_module_cre = [node for node in flow_network.successors(module_cre) if node in all_outflow_vars]

                    if len(outflow_successors_of_module_cre) != 0 and module_cre not in module_cre_vars:
                        module_cre_vars.append(module_cre)
                            
                        for outflow in outflow_successor_of_module_cre:
                            if outflow not in outflow_vars:
                                outflow_vars.append(outflow)
                
            # If outflow vars have been specified, we get the GEMs connected to the specified outflows, and then the inflow connected to the inferred GEMs
            else:
                inflow_vars = []
                module_tf_vars = []
                module_cre_vars = []
                
                for outflow in outflow_vars:
                    outflow_predecessors = [node for node in flow_network.predecessors(outflow) if node in all_module_cre_vars]

                    for module_cre in outflow_predecessors:
                        module_cre_predecessors = [node for node in flow_network.predecessors(module_cre) if node in all_module_tf_vars]
                        if len(module_cre_predecessors) != 0 and module_cre not in module_cre_vars:
                            module_cre_vars.append(module_cre)

                

                for module_cre in module_cre_vars:
                    module_cre_predecessors = [node for node in flow_network.predecessors(module_cre) if node in all_module_tf_vars]

                    for module_tf in module_cre_predecessors:
                        if module_tf not in module_tf_vars:
                            module_tf_vars.append(module_tf)

    elif var_types_used.sum() == 2:

        # We need the downstream CRE modules that are downstream of a TF module and connected to an outflow node as well
        if inflows_used and module_tfs_used:
            module_cre_vars = list(set([edge[1] for edge in flow_network.edges() if edge[0] in module_tf_vars and edge[1] in all_module_cre_vars\
                                         and len([succ for succ in flow_network.successors(edge[1]) if succ in all_outflow_vars]) > 0]))
            outflow_vars = list(set([edge[1] for edge in flow_network.edges() if edge[0] in module_cre_vars and edge[1] in all_outflow_vars]))

        # We need the module TF variables that connect the inflows to the CREs, and the downstream outflow nodfes of the CRE modules
        elif inflows_used and module_cres_used:
            module_tf_vars = list(set([edge[0] for edge in flow_network.edges() if edge[1] in module_cre_vars and edge[0] in all_module_tf_vars\
                                            and len([pred for pred in flow_network.predecessors(edge[0]) if pred in all_inflow_vars]) > 0]))
            outflow_vars = list(set([edge[1] for edge in flow_network.edges() if edge[0] in module_cre_vars and edge[1] in all_outflow_vars]))
        
        # We need the module TFs downstream of the inflows that are connected to CRE modules that are upstream of the outflow_nodes
        
        elif inflows_used and outflows_used:
            module_tfs_from_inflows = list(set([edge[1] for edge in flow_network.edges() if edge[0] in inflow_vars and edge[1] in all_module_tf_vars]))
            module_cres_to_outflows = list(set([edge[0] for edge in flow_network.edges() if edge[0] in all_module_cre_vars and edge[1] in outflow_vars]))
            
            for module_tf in module_tfs_from_inflows:
                for module_cre in module_cres_to_outflows:
                    if flow_network.has_edge(module_tf, module_cre):
                        if module_tf not in module_tf_vars:
                            module_tf_vars.append(module_tf)
                        if module_cre not in module_cre_vars:
                            module_cre_vars.append(module_cre)

        # Get the inflow nodes connecting to a TF module and the outflow nodes downstream of a CRE module
        elif module_tfs_used and module_cres_used:
            inflow_vars = list(set([edge[0] for edge in flow_network.edges() if edge[1] in module_tf_vars and edge[0] in all_inflow_vars]))
            outflow_vars = list(set([edge[1] for edge in flow_network.edges() if edge[0] in module_cre_vars and edge[1] in all_outflow_vars]))

        # Get the CRE modules that connect a specified TF module to a specified outflow node
        elif module_tfs_used and outflows_used:
            inflow_vars = list(set([edge[0] for edge in flow_network.edges() if edge[1] in module_tf_vars and edge[0] in all_inflow_vars]))

            module_cres_from_module_tfs = list(set([edge[1] for edge in flow_network.edges() if edge[0] in module_tf_vars and edge[1] in all_module_cre_vars]))
            module_cres_to_outflows = list(set([edge[0] for edge in flow_network.edges() if edge[0] in all_module_cre_vars and edge[1] in outflow_vars]))
            module_cre_vars = list(set(module_cres_from_module_tfs) & set(module_cres_to_outflows)) # We want the intersection

        # Get the TF modules connected to the CRE modules, provided the TF module has an upstream inflow node connecting, 
        # and then the inflows connected to those TF modules
        elif module_cres_used and outflows_used:
            module_tf_vars = list(set([edge[0] for edge in flow_network.edges() if edge[1] in module_cre_vars and edge[0] in all_module_tf_vars\
                                            and len([pred for pred in flow_network.predecessors(edge[0]) if pred in all_inflow_vars]) > 0]))
            inflow_vars = list(set([edge[0] for edge in flow_network.edges() if edge[1] in module_tf_vars and edge[0] in all_inflow_vars]))

    elif var_types_used.sum() == 3:
        
        if not inflows_used:
            inflow_vars = list(set([edge[0] for edge in flow_network.edges() if edge[0] in all_inflow_vars and edge[1] in module_tf_vars]))

        elif not module_tfs_used:
            module_tf_vars = list(set([edge[1] for edge in flow_network.edges() if edge[0] in inflow_vars and edge[1] in all_module_tf_vars and edge[1] in module_cre_vars]))

        elif not module_cres_used:
            module_cre_vars = list(set([edge[1] for edge in flow_network.edges() if edge[0] in module_tf_vars and edge[1] in all_module_cre_vars and edge[1] in outflow_vars]))

        else: # Not outflows used
            outflow_vars = list(set([edge[1] for edge in flow_network.edges() if edge[0] in module_cre_vars and edge[1] in all_outflow_vars]))

    # We will sort the inflow vars etc now
    inflow_vars = sorted(inflow_vars)
    module_tf_vars = sorted(module_tf_vars, key=lambda x: int(x.split('-')[1])) # Assumption is that all module vars are written like 'GEM-1', etc.
    module_cre_vars = sorted(module_cre_vars, key=lambda x: int(x.split('-')[1])) # Assumption is that all module vars are written like 'GEM-1', etc.

    outflow_vars = sorted(outflow_vars)
    resultant_pattern_graph = flow_network.subgraph(inflow_vars + module_tf_vars + module_cre_vars + outflow_vars)

    node_colours = {}
    count = 0

    for inflow_var in all_inflow_vars:
        node_colours[inflow_var] = palette_network[count]
        count += 1

    for module_var in all_module_tf_vars:
        node_colours[module_var] = palette_network[count]
        count += 1
    for module_var in all_module_cre_vars:
        node_colours[module_var] = palette_network[count]
        count += 1
        
    for ouflow_var in all_outflow_vars:
        node_colours[ouflow_var] = palette_network[count]
        count += 1

    for node in resultant_pattern_graph.nodes():

        if align_mode == 'horizontal':

            if node in inflow_vars:
                resultant_pattern_graph.nodes[node]['level'] = 4

            elif node in module_tf_vars:
                resultant_pattern_graph.nodes[node]['level'] = 3

            elif node in module_cre_vars:
                resultant_pattern_graph.nodes[node]['level'] = 2

            else: # Downstream outflow vars
                resultant_pattern_graph.nodes[node]['level'] = 1

        else: # Otherwise we're going veritcally

            if node in inflow_vars:
                resultant_pattern_graph.nodes[node]['level'] = 1

            elif node in module_tf_vars:
                resultant_pattern_graph.nodes[node]['level'] = 2
            
            elif node in module_cre_vars:
                resultant_pattern_graph.nodes[node]['level'] = 3

            else:
                resultant_pattern_graph.nodes[node]['level'] = 4

    # We will sort everything one last time out of paranoia
    resultant_pattern_graph_pos = nx.multipartite_layout(resultant_pattern_graph, subset_key='level', align=align_mode, scale=2.0)

    # Fix the damned positions
    inflow_pos = [resultant_pattern_graph_pos[node] for node in inflow_vars]
    module_tf_pos = [resultant_pattern_graph_pos[node] for node in module_tf_vars]
    module_cre_pos = [resultant_pattern_graph_pos[node] for node in module_cre_vars]
    outflow_pos = [resultant_pattern_graph_pos[node] for node in outflow_vars]

    if align_mode == 'horizontal':
        inflow_pos = sorted(inflow_pos, key=lambda x: x[0], reverse=False)
        module_tf_pos = sorted(module_tf_pos, key=lambda x: x[0], reverse=False)
        module_cre_pos = sorted(module_cre_pos, key=lambda x: x[0], reverse=False)
        outflow_pos = sorted(outflow_pos, key=lambda x: x[0], reverse=False)

    else:
        inflow_pos = sorted(inflow_pos, key=lambda x: x[1], reverse=False)
        module_tf_pos = sorted(module_tf_pos, key=lambda x: x[1], reverse=False)
        module_cre_pos = sorted(module_cre_pos, key=lambda x: x[1], reverse=False)
        outflow_pos = sorted(outflow_pos, key=lambda x: x[1], reverse=False)

    for i, node in enumerate(inflow_vars):
        resultant_pattern_graph_pos[node] = inflow_pos[i]
        
    for i, node in enumerate(module_tf_vars):
        resultant_pattern_graph_pos[node] = module_tf_pos[i]
        
    for i, node in enumerate(module_cre_vars):
        resultant_pattern_graph_pos[node] = module_cre_pos[i]

    for i, node in enumerate(outflow_vars):
        resultant_pattern_graph_pos[node] = outflow_pos[i]

    # Let's sort the nodes by alphabetical order
    resultant_pattern_graph_edge_colours = [node_colours[edge[0]] for edge in resultant_pattern_graph.edges()]
    edge_widths = [width_scale*resultant_pattern_graph[u][v]['weight'] for u,v in resultant_pattern_graph.edges()]

    if ax is None:
        fig, ax = plt.subplots()
    nx.draw_networkx_edges(resultant_pattern_graph,
                            resultant_pattern_graph_pos,
                            edge_color=resultant_pattern_graph_edge_colours,
                            width=edge_widths,
                             alpha=0.75,
                             connectionstyle="arc3,rad=0.2",
                             ax=ax)
    nx.draw_networkx_nodes(resultant_pattern_graph,
                           resultant_pattern_graph_pos,
                           nodelist=inflow_vars, node_color=[node_colours[node] for node in inflow_vars], linewidths=1.0, edgecolors='black', node_size=100, ax=ax)
    nx.draw_networkx_nodes(resultant_pattern_graph,
                           resultant_pattern_graph_pos,
                           nodelist=module_tf_vars, node_color=[node_colours[node] for node in module_tf_vars], node_shape='v', node_size=100, ax=ax)
    nx.draw_networkx_nodes(resultant_pattern_graph,
                           resultant_pattern_graph_pos,
                           nodelist=module_cre_vars, node_color=[node_colours[node] for node in module_cre_vars], node_shape='d', node_size=100, ax=ax)
    nx.draw_networkx_nodes(resultant_pattern_graph,
                           resultant_pattern_graph_pos,
                           nodelist=outflow_vars, node_color=[node_colours[node] for node in outflow_vars], node_shape='s', node_size=100, ax=ax)

    inflow_labels = nx.draw_networkx_labels(resultant_pattern_graph,
                                            resultant_pattern_graph_pos,
                                            labels = {node:node for node in resultant_pattern_graph.nodes() if node in inflow_vars},
                                            font_size=14,
                                            font_family='Arial',
                                            ax=ax,
                                            horizontalalignment='left',
                                            verticalalignment = 'bottom');
    for text in inflow_labels.values():
        text.set_rotation(45)

    module_tf_labels = nx.draw_networkx_labels(resultant_pattern_graph,
                                            resultant_pattern_graph_pos,
                                            labels = {node:node for node in resultant_pattern_graph.nodes() if node in module_tf_vars},
                                            font_size=14, font_family='Arial',
                                            ax=ax,
                                            horizontalalignment='left',
                                            verticalalignment = 'bottom');
    for text in module_tf_labels.values():
        text.set_rotation(45)

    module_cre_labels = nx.draw_networkx_labels(resultant_pattern_graph,
                                            resultant_pattern_graph_pos,
                                            labels = {node:node for node in resultant_pattern_graph.nodes() if node in module_cre_vars},
                                            font_size=14, font_family='Arial',
                                            ax=ax,
                                            horizontalalignment='left',
                                            verticalalignment = 'bottom');
    for text in module_cre_labels.values():
        text.set_rotation(45)
        
    outflow_labels = nx.draw_networkx_labels(resultant_pattern_graph,
                                            resultant_pattern_graph_pos,
                                            labels = {node:node for node in resultant_pattern_graph.nodes() if node in outflow_vars},
                                            font_size=14, font_family='Arial',
                                            ax=ax,
                                            horizontalalignment='right',
                                            verticalalignment = 'top');
    for text in outflow_labels.values():
        text.set_rotation(45)

    plt.tight_layout()
    plt.margins(x=x_margin_offset)
    plt.margins(y=y_margin_offset)
    plt.box(False)

    return ax