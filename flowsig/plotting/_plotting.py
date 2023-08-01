import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import scanpy as sc
from typing import Union, Sequence

palette_network = list(sns.color_palette("tab20") \
                       + sns.color_palette("tab20b") \
                       + sns.color_palette("tab20c")\
                       + sns.color_palette("Set1")\
                       + sns.color_palette("Set2"))

# def plot_differentially_inflowing_signals():

# def plot_differentially_outflowing_signals():

def plot_intercellular_flows_from_inflows(adata: sc.AnnData,
                                          inflow_vars: Union[str, Sequence[str]],
                                          flow_network: nx.DiGraph(),
                                          flowsig_network_key: str = 'flowsig_network',
                                          align_mode: str = 'horizontal',
                                          width_scale: float = 2.0,
                                          x_margin_offset: float = 0.3,
                                          y_margin_offset: float = 0.0):
    
    flow_vars = adata.uns[flowsig_network_key]['flow_var_info'].index.tolist()
    flow_var_info = adata.uns[flowsig_network_key]['flow_var_info']

    all_gems = [node for node in flow_vars\
                        if flow_var_info.loc[node]['Type'] == 'module']
    all_outflow_vars = [node for node in flow_vars\
                        if flow_var_info.loc[node]['Type'] == 'outflow']

    node_colours = {}
    count = 0

    for inflow_var in inflow_vars:
        node_colours[inflow_var] = palette_network[count]
        count += 1

    associated_gems = sorted([edge[1] for edge in flow_network.edges() if edge[0] in inflow_vars and edge[1] in all_gems])
    for gem in associated_gems:
        node_colours[gem] = palette_network[count]
        count += 1

    downstream_outflows = []
    for gem in associated_gems:

        outflow_successors = [node for node in flow_network.successors(gem) if node in all_outflow_vars]
        
        for outflow_var in outflow_successors:
            if outflow_var not in downstream_outflows:
                downstream_outflows.append(outflow_var)            

    downstream_outflows = sorted(downstream_outflows)

    for outflow_var in downstream_outflows:
        node_colours[outflow_var] = palette_network[count]
        count += 1

    print(inflow_vars + associated_gems + downstream_outflows)
    resultant_pattern_graph = flow_network.subgraph(inflow_vars + associated_gems + downstream_outflows)

    for node in resultant_pattern_graph.nodes():

        if align_mode == 'horizontal':

            if node in inflow_vars:
                resultant_pattern_graph.nodes[node]['level'] = 3

            elif node in associated_gems:
                resultant_pattern_graph.nodes[node]['level'] = 2

            else: # Downstream outflow vars
                resultant_pattern_graph.nodes[node]['level'] = 1

        else: # Otherwise we're going veritcally

            if node in inflow_vars:
                resultant_pattern_graph.nodes[node]['level'] = 1

            elif node in associated_gems:
                resultant_pattern_graph.nodes[node]['level'] = 2

            else:
                resultant_pattern_graph.nodes[node]['level'] = 3

    resultant_pattern_graph_pos = nx.multipartite_layout(resultant_pattern_graph, subset_key='level', align=align_mode, scale=2.0)

    print(resultant_pattern_graph_pos)

    resultant_pattern_graph_edge_colours = [node_colours[edge[0]] for edge in resultant_pattern_graph.edges()]
    edge_widths = [width_scale*resultant_pattern_graph[u][v]['weight'] for u,v in resultant_pattern_graph.edges()]
    nx.draw_networkx_edges(resultant_pattern_graph, resultant_pattern_graph_pos, edge_color=resultant_pattern_graph_edge_colours, width=edge_widths, alpha=0.75, connectionstyle="arc3,rad=0.2")
    nx.draw_networkx_labels(resultant_pattern_graph, resultant_pattern_graph_pos, font_size=14, font_family='Arial');
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=inflow_vars, node_color=[node_colours[node] for node in inflow_vars], linewidths=1.0, edgecolors='black', node_size=500)
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=associated_gems, node_color=[node_colours[node] for node in associated_gems], node_shape='v', node_size=100)
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=downstream_outflows, node_color=[node_colours[node] for node in downstream_outflows], node_shape='s', node_size=100)

    plt.tight_layout()
    plt.margins(x=x_margin_offset)
    plt.margins(y=y_margin_offset)
    plt.box(False)

def plot_intercellular_flows_from_gems(adata: sc.AnnData,
                                          gem_vars: Union[str, Sequence[str]],
                                          flow_network: nx.DiGraph(),
                                          flowsig_network_key: str = 'flowsig_network',
                                          align_mode: str = 'horizontal',
                                          width_scale: float = 2.0,
                                          x_margin_offset: float = 0.3,
                                          y_margin_offset: float = 0.0):
    
    flow_vars = adata.uns[flowsig_network_key]['flow_var_info'].index.tolist()
    flow_var_info = adata.uns[flowsig_network_key]['flow_var_info']

    all_inflow_vars = [node for node in flow_vars\
                        if flow_var_info.loc[node]['Type'] == 'inflow']
    all_outflow_vars = [node for node in flow_vars\
                        if flow_var_info.loc[node]['Type'] == 'outflow']

    node_colours = {}
    count = 0

    upstream_inflows = []
    for gem in gem_vars:

        inflow_predecessors = [node for node in flow_network.predecessors(gem) if node in all_inflow_vars]
        
        for inflow_var in inflow_predecessors:
            if inflow_var not in upstream_inflows:
                upstream_inflows.append(inflow_var)            

    upstream_inflows = sorted(upstream_inflows)

    for inflow_var in upstream_inflows:
        node_colours[inflow_var] = palette_network[count]
        count += 1

    downstream_outflows = []
    for gem in gem_vars:

        outflow_successors = [node for node in flow_network.successors(gem) if node in all_outflow_vars]
        
        for outflow_var in outflow_successors:
            if outflow_var not in downstream_outflows:
                downstream_outflows.append(outflow_var)            

    downstream_outflows = sorted(downstream_outflows)

    for outflow_var in downstream_outflows:
        node_colours[outflow_var] = palette_network[count]
        count += 1

    resultant_pattern_graph = flow_network.subgraph(upstream_inflows + gem_vars + downstream_outflows)

    for node in resultant_pattern_graph.nodes():

        if align_mode == 'horizontal':

            if node in upstream_inflows:
                resultant_pattern_graph.nodes[node]['level'] = 3

            elif node in gem_vars:
                resultant_pattern_graph.nodes[node]['level'] = 2

            else: # Downstream outflow vars
                resultant_pattern_graph.nodes[node]['level'] = 1

        else: # Otherwise we're going veritcally

            if node in upstream_inflows:
                resultant_pattern_graph.nodes[node]['level'] = 1

            elif node in gem_vars:
                resultant_pattern_graph.nodes[node]['level'] = 2

            else:
                resultant_pattern_graph.nodes[node]['level'] = 3

    resultant_pattern_graph_pos = nx.multipartite_layout(resultant_pattern_graph, subset_key='level', align=align_mode, scale=2.0)

    resultant_pattern_graph_edge_colours = [node_colours[edge[0]] for edge in resultant_pattern_graph.edges()]
    edge_widths = [width_scale*resultant_pattern_graph[u][v]['weight'] for u,v in resultant_pattern_graph.edges()]
    nx.draw_networkx_edges(resultant_pattern_graph, resultant_pattern_graph_pos, edge_color=resultant_pattern_graph_edge_colours, width=edge_widths, alpha=0.75, connectionstyle="arc3,rad=0.2")
    nx.draw_networkx_labels(resultant_pattern_graph, resultant_pattern_graph_pos, font_size=14, font_family='Arial');
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=upstream_inflows, node_color=[node_colours[node] for node in upstream_inflows], linewidths=1.0, edgecolors='black', node_size=100)
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=gem_vars, node_color=[node_colours[node] for node in gem_vars], node_shape='v', node_size=500)
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=downstream_outflows, node_color=[node_colours[node] for node in downstream_outflows], node_shape='s', node_size=100)

    plt.tight_layout()
    plt.margins(x=x_margin_offset)
    plt.margins(y=y_margin_offset)
    plt.box(False)

def plot_intercellular_flows_from_outflows(adata: sc.AnnData,
                                          outflow_vars: Union[str, Sequence[str]],
                                          flow_network: nx.DiGraph(),
                                          flowsig_network_key: str = 'flowsig_network',
                                          align_mode: str = 'horizontal',
                                          width_scale: float = 2.0,
                                          x_margin_offset: float = 0.3,
                                          y_margin_offset: float = 0.0):
    
    flow_vars = adata.uns[flowsig_network_key]['flow_var_info'].index.tolist()
    flow_var_info = adata.uns[flowsig_network_key]['flow_var_info']

    all_inflow_vars = [node for node in flow_vars\
                        if flow_var_info.loc[node]['Type'] == 'inflow']
    all_gems = [node for node in flow_vars\
                        if flow_var_info.loc[node]['Type'] == 'module']

    node_colours = {}
    count = 0

    for outflow_var in outflow_vars:
        node_colours[outflow_var] = palette_network[count]
        count += 1

    associated_gems = sorted([edge[0] for edge in flow_network.edges() if edge[1] in outflow_vars and edge[0] in all_gems])
    for gem in associated_gems:
        node_colours[gem] = palette_network[count]
        count += 1

    upstream_inflows = []
    for gem in associated_gems:

        inflow_predecessors = [node for node in flow_network.predecessors(gem) if node in all_inflow_vars]
        
        for inflow_var in inflow_predecessors:
            if inflow_var not in upstream_inflows:
                upstream_inflows.append(inflow_var)            

    upstream_inflows = sorted(upstream_inflows)

    for inflow_var in upstream_inflows:
        node_colours[outflow_var] = palette_network[count]
        count += 1

    resultant_pattern_graph = flow_network.subgraph(upstream_inflows + associated_gems + outflow_vars)

    for node in resultant_pattern_graph.nodes():

        if align_mode == 'horizontal':

            if node in upstream_inflows:
                resultant_pattern_graph.nodes[node]['level'] = 3

            elif node in associated_gems:
                resultant_pattern_graph.nodes[node]['level'] = 2

            else: # Downstream outflow vars
                resultant_pattern_graph.nodes[node]['level'] = 1

        else: # Otherwise we're going veritcally

            if node in upstream_inflows:
                resultant_pattern_graph.nodes[node]['level'] = 1

            elif node in associated_gems:
                resultant_pattern_graph.nodes[node]['level'] = 2

            else:
                resultant_pattern_graph.nodes[node]['level'] = 3

    resultant_pattern_graph_pos = nx.multipartite_layout(resultant_pattern_graph, subset_key='level', align=align_mode, scale=2.0)

    resultant_pattern_graph_edge_colours = [node_colours[edge[0]] for edge in resultant_pattern_graph.edges()]
    edge_widths = [width_scale*resultant_pattern_graph[u][v]['weight'] for u,v in resultant_pattern_graph.edges()]
    nx.draw_networkx_edges(resultant_pattern_graph, resultant_pattern_graph_pos, edge_color=resultant_pattern_graph_edge_colours, width=edge_widths, alpha=0.75, connectionstyle="arc3,rad=0.2")
    nx.draw_networkx_labels(resultant_pattern_graph, resultant_pattern_graph_pos, font_size=14, font_family='Arial');
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=upstream_inflows, node_color=[node_colours[node] for node in upstream_inflows], linewidths=1.0, edgecolors='black', node_size=100)
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=associated_gems, node_color=[node_colours[node] for node in associated_gems], node_shape='v', node_size=100)
    nx.draw_networkx_nodes(resultant_pattern_graph, resultant_pattern_graph_pos, nodelist=outflow_vars, node_color=[node_colours[node] for node in outflow_vars], node_shape='s', node_size=500)

    plt.tight_layout()
    plt.margins(x=x_margin_offset)
    plt.margins(y=y_margin_offset)
    plt.box(False)
