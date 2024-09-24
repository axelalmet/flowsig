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
                       + sns.color_palette("Set2"))

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

        ax.text(point['x']+0.3, point['y'] + 0.1, str(point['val']), fontdict={'size':12.0})


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
    result["-logQ"] = -np.log(result["pvals_adj"].astype("float"))

    lowqval_de = result.loc[(abs(result["logfoldchanges"]) > logfc_threshold)&(abs(result["pvals_adj"]) < qval_threshold)]
    other_de = result.loc[(abs(result["logfoldchanges"]) <= logfc_threshold)|(abs(result["pvals_adj"]) >= qval_threshold)]

    if ax is None:
        fig, ax = plt.subplots()

    sns.regplot(
        x=other_de["logfoldchanges"],
        y=other_de["-logQ"],
        fit_reg=False,
        scatter_kws={"s": scatter_size},
        ax=ax
    )
    sns.regplot(
        x=lowqval_de["logfoldchanges"],
        y=lowqval_de["-logQ"],
        fit_reg=False,
        scatter_kws={"s": scatter_size},
        ax=ax
    )
    ax.set_xlabel("log2 FC")
    ax.set_ylabel("-log Q-value")

    if label_lowqval:
        label_point(lowqval_de['logfoldchanges'],
                    lowqval_de['-logQ'],
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

    elif var_types_used.sum(2) == 2:

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
    inflow_vars = sorted(inflow_vars)
    module_vars = sorted(module_vars, key=lambda x: int(x[4:])) # Assumption is that all module vars are written like 'GEM-1', etc.
    outflow_vars = sorted(outflow_vars)
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
    nx.draw_networkx_edges(resultant_pattern_graph, resultant_pattern_graph_pos, edge_color=resultant_pattern_graph_edge_colours, width=edge_widths, alpha=0.75, connectionstyle="arc3,rad=0.2", ax=ax)
    nx.draw_networkx_labels(resultant_pattern_graph, resultant_pattern_graph_pos, font_size=14, font_family='Arial', ax=ax);
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=inflow_vars, node_color=[node_colours[node] for node in inflow_vars], linewidths=1.0, edgecolors='black', node_size=500, ax=ax)
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=module_vars, node_color=[node_colours[node] for node in module_vars], node_shape='v', node_size=100, ax=ax)
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=outflow_vars, node_color=[node_colours[node] for node in outflow_vars], node_shape='s', node_size=100, ax=ax)

    plt.tight_layout()
    plt.margins(x=x_margin_offset)
    plt.margins(y=y_margin_offset)
    plt.box(False)

    return ax

# def plot_intercellular_flows_from_inflows(adata: sc.AnnData,
#                                           inflow_vars: Union[str, Sequence[str]],
#                                           flow_network: nx.DiGraph,
#                                           flowsig_network_key: str = 'flowsig_network',
#                                           align_mode: str = 'horizontal',
#                                           width_scale: float = 2.0,
#                                           x_margin_offset: float = 0.3,
#                                           y_margin_offset: float = 0.0):
    
#     flow_vars = adata.uns[flowsig_network_key]['flow_var_info'].index.tolist()
#     flow_var_info = adata.uns[flowsig_network_key]['flow_var_info']

#     # inflow_vars = []
#     all_gems = [node for node in flow_vars\
#                         if flow_var_info.loc[node]['Type'] == 'module']
#     all_outflow_vars = [node for node in flow_vars\
#                         if flow_var_info.loc[node]['Type'] == 'outflow']

#     node_colours = {}
#     count = 0

#     for inflow_var in inflow_vars:
#         node_colours[inflow_var] = palette_network[count]
#         count += 1

#     associated_gems = sorted([edge[1] for edge in flow_network.edges() if edge[0] in inflow_vars and edge[1] in all_gems])
#     for gem in associated_gems:
#         node_colours[gem] = palette_network[count]
#         count += 1

#     downstream_outflows = []
#     for gem in associated_gems:

#         outflow_successors = [node for node in flow_network.successors(gem) if node in all_outflow_vars]
        
#         for outflow_var in outflow_successors:
#             if outflow_var not in downstream_outflows:
#                 downstream_outflows.append(outflow_var)            

#     downstream_outflows = sorted(downstream_outflows)

#     for outflow_var in downstream_outflows:
#         node_colours[outflow_var] = palette_network[count]
#         count += 1

#     resultant_pattern_graph = flow_network.subgraph(inflow_vars + associated_gems + downstream_outflows)

#     for node in resultant_pattern_graph.nodes():

#         if align_mode == 'horizontal':

#             if node in inflow_vars:
#                 resultant_pattern_graph.nodes[node]['level'] = 3

#             elif node in associated_gems:
#                 resultant_pattern_graph.nodes[node]['level'] = 2

#             else: # Downstream outflow vars
#                 resultant_pattern_graph.nodes[node]['level'] = 1

#         else: # Otherwise we're going veritcally

#             if node in inflow_vars:
#                 resultant_pattern_graph.nodes[node]['level'] = 1

#             elif node in associated_gems:
#                 resultant_pattern_graph.nodes[node]['level'] = 2

#             else:
#                 resultant_pattern_graph.nodes[node]['level'] = 3

#     resultant_pattern_graph_pos = nx.multipartite_layout(resultant_pattern_graph, subset_key='level', align=align_mode, scale=2.0)

#     # Let's sort the nodes by alphabetical order
#     resultant_pattern_graph_edge_colours = [node_colours[edge[0]] for edge in resultant_pattern_graph.edges()]
#     edge_widths = [width_scale*resultant_pattern_graph[u][v]['weight'] for u,v in resultant_pattern_graph.edges()]

#     fig, ax = plt.subplots()
#     nx.draw_networkx_edges(resultant_pattern_graph, resultant_pattern_graph_pos, edge_color=resultant_pattern_graph_edge_colours, width=edge_widths, alpha=0.75, connectionstyle="arc3,rad=0.2", ax=ax)
#     nx.draw_networkx_labels(resultant_pattern_graph, resultant_pattern_graph_pos, font_size=14, font_family='Arial', ax=ax);
#     nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=inflow_vars, node_color=[node_colours[node] for node in inflow_vars], linewidths=1.0, edgecolors='black', node_size=500, ax=ax)
#     nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=associated_gems, node_color=[node_colours[node] for node in associated_gems], node_shape='v', node_size=100, ax=ax)
#     nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=downstream_outflows, node_color=[node_colours[node] for node in downstream_outflows], node_shape='s', node_size=100, ax=ax)

#     plt.tight_layout()
#     plt.margins(x=x_margin_offset)
#     plt.margins(y=y_margin_offset)
#     plt.box(False)

#     return ax

# def plot_intercellular_flows_from_gems(adata: sc.AnnData,
#                                           gem_vars: Union[str, Sequence[str]],
#                                           flow_network: nx.DiGraph,
#                                           flowsig_network_key: str = 'flowsig_network',
#                                           align_mode: str = 'horizontal',
#                                           width_scale: float = 2.0,
#                                           x_margin_offset: float = 0.3,
#                                           y_margin_offset: float = 0.0):
    
#     flow_vars = adata.uns[flowsig_network_key]['flow_var_info'].index.tolist()
#     flow_var_info = adata.uns[flowsig_network_key]['flow_var_info']

#     all_inflow_vars = [node for node in flow_vars\
#                         if flow_var_info.loc[node]['Type'] == 'inflow']
#     all_outflow_vars = [node for node in flow_vars\
#                         if flow_var_info.loc[node]['Type'] == 'outflow']

#     node_colours = {}
#     count = 0

#     upstream_inflows = []
#     for gem in gem_vars:

#         inflow_predecessors = [node for node in flow_network.predecessors(gem) if node in all_inflow_vars]
        
#         for inflow_var in inflow_predecessors:
#             if inflow_var not in upstream_inflows:
#                 upstream_inflows.append(inflow_var)            

#     upstream_inflows = sorted(upstream_inflows)

#     for inflow_var in upstream_inflows:
#         node_colours[inflow_var] = palette_network[count]
#         count += 1

#     downstream_outflows = []
#     for gem in gem_vars:

#         outflow_successors = [node for node in flow_network.successors(gem) if node in all_outflow_vars]
        
#         for outflow_var in outflow_successors:
#             if outflow_var not in downstream_outflows:
#                 downstream_outflows.append(outflow_var)            

#     downstream_outflows = sorted(downstream_outflows)

#     for outflow_var in downstream_outflows:
#         node_colours[outflow_var] = palette_network[count]
#         count += 1

#     resultant_pattern_graph = flow_network.subgraph(upstream_inflows + gem_vars + downstream_outflows)

#     for node in resultant_pattern_graph.nodes():

#         if align_mode == 'horizontal':

#             if node in upstream_inflows:
#                 resultant_pattern_graph.nodes[node]['level'] = 3

#             elif node in gem_vars:
#                 resultant_pattern_graph.nodes[node]['level'] = 2

#             else: # Downstream outflow vars
#                 resultant_pattern_graph.nodes[node]['level'] = 1

#         else: # Otherwise we're going veritcally

#             if node in upstream_inflows:
#                 resultant_pattern_graph.nodes[node]['level'] = 1

#             elif node in gem_vars:
#                 resultant_pattern_graph.nodes[node]['level'] = 2

#             else:
#                 resultant_pattern_graph.nodes[node]['level'] = 3

#     resultant_pattern_graph_pos = nx.multipartite_layout(resultant_pattern_graph, subset_key='level', align=align_mode, scale=2.0)

#     resultant_pattern_graph_edge_colours = [node_colours[edge[0]] for edge in resultant_pattern_graph.edges()]
#     edge_widths = [width_scale*resultant_pattern_graph[u][v]['weight'] for u,v in resultant_pattern_graph.edges()]
#     nx.draw_networkx_edges(resultant_pattern_graph, resultant_pattern_graph_pos, edge_color=resultant_pattern_graph_edge_colours, width=edge_widths, alpha=0.75, connectionstyle="arc3,rad=0.2")
#     nx.draw_networkx_labels(resultant_pattern_graph, resultant_pattern_graph_pos, font_size=14, font_family='Arial');
#     nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=upstream_inflows, node_color=[node_colours[node] for node in upstream_inflows], linewidths=1.0, edgecolors='black', node_size=100)
#     nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=gem_vars, node_color=[node_colours[node] for node in gem_vars], node_shape='v', node_size=500)
#     nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=downstream_outflows, node_color=[node_colours[node] for node in downstream_outflows], node_shape='s', node_size=100)

#     plt.tight_layout()
#     plt.margins(x=x_margin_offset)
#     plt.margins(y=y_margin_offset)
#     plt.box(False)
#     plt.show()

# def plot_intercellular_flows_from_outflows(adata: sc.AnnData,
#                                           outflow_vars: Union[str, Sequence[str]],
#                                           flow_network: nx.DiGraph,
#                                           flowsig_network_key: str = 'flowsig_network',
#                                           align_mode: str = 'horizontal',
#                                           width_scale: float = 2.0,
#                                           x_margin_offset: float = 0.3,
#                                           y_margin_offset: float = 0.0):
    
#     flow_vars = adata.uns[flowsig_network_key]['flow_var_info'].index.tolist()
#     flow_var_info = adata.uns[flowsig_network_key]['flow_var_info']

#     all_inflow_vars = [node for node in flow_vars\
#                         if flow_var_info.loc[node]['Type'] == 'inflow']
#     all_gems = [node for node in flow_vars\
#                         if flow_var_info.loc[node]['Type'] == 'module']

#     node_colours = {}
#     count = 0

#     for outflow_var in outflow_vars:
#         node_colours[outflow_var] = palette_network[count]
#         count += 1

#     associated_gems = sorted([edge[0] for edge in flow_network.edges() if edge[1] in outflow_vars and edge[0] in all_gems])
#     for gem in associated_gems:
#         node_colours[gem] = palette_network[count]
#         count += 1

#     upstream_inflows = []
#     for gem in associated_gems:

#         inflow_predecessors = [node for node in flow_network.predecessors(gem) if node in all_inflow_vars]
        
#         for inflow_var in inflow_predecessors:
#             if inflow_var not in upstream_inflows:
#                 upstream_inflows.append(inflow_var)            

#     upstream_inflows = sorted(upstream_inflows)

#     for inflow_var in upstream_inflows:
#         node_colours[outflow_var] = palette_network[count]
#         count += 1

#     resultant_pattern_graph = flow_network.subgraph(upstream_inflows + associated_gems + outflow_vars)

#     for node in resultant_pattern_graph.nodes():

#         if align_mode == 'horizontal':

#             if node in upstream_inflows:
#                 resultant_pattern_graph.nodes[node]['level'] = 3

#             elif node in associated_gems:
#                 resultant_pattern_graph.nodes[node]['level'] = 2

#             else: # Downstream outflow vars
#                 resultant_pattern_graph.nodes[node]['level'] = 1

#         else: # Otherwise we're going veritcally

#             if node in upstream_inflows:
#                 resultant_pattern_graph.nodes[node]['level'] = 1

#             elif node in associated_gems:
#                 resultant_pattern_graph.nodes[node]['level'] = 2

#             else:
#                 resultant_pattern_graph.nodes[node]['level'] = 3

#     resultant_pattern_graph_pos = nx.multipartite_layout(resultant_pattern_graph, subset_key='level', align=align_mode, scale=2.0)

#     resultant_pattern_graph_edge_colours = [node_colours[edge[0]] for edge in resultant_pattern_graph.edges()]
#     edge_widths = [width_scale*resultant_pattern_graph[u][v]['weight'] for u,v in resultant_pattern_graph.edges()]
#     nx.draw_networkx_edges(resultant_pattern_graph, resultant_pattern_graph_pos, edge_color=resultant_pattern_graph_edge_colours, width=edge_widths, alpha=0.75, connectionstyle="arc3,rad=0.2")
#     nx.draw_networkx_labels(resultant_pattern_graph, resultant_pattern_graph_pos, font_size=14, font_family='Arial');
#     nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=upstream_inflows, node_color=[node_colours[node] for node in upstream_inflows], linewidths=1.0, edgecolors='black', node_size=100)
#     nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=associated_gems, node_color=[node_colours[node] for node in associated_gems], node_shape='v', node_size=100)
#     nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=outflow_vars, node_color=[node_colours[node] for node in outflow_vars], node_shape='s', node_size=500)

#     plt.tight_layout()
#     plt.margins(x=x_margin_offset)
#     plt.margins(y=y_margin_offset)
#     plt.box(False)
#     plt.show()
