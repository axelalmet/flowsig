from typing import Optional, Callable, Tuple
from anndata import AnnData
import networkx as nx
import numpy as np
import pandas as pd
import graphical_models as gpm

def _total_edge_weights(
        adjacency: np.ndarray
) -> dict[Tuple[int, int], float]:
    nz_rows, nz_cols = adjacency.nonzero()
    edge_weights: dict[Tuple[int, int], float] = {}
    for i, j in zip(nz_rows, nz_cols):
        a, b = (i, j) if i < j else (j, i)
        edge_weights[(a, b)] = edge_weights.get((a, b), 0.0) + adjacency[i, j]
    return edge_weights

def _filter_edges_by_weight(
    cpdag: gpm.PDAG,
    adjacency: np.ndarray,
    keep_if: Callable[[int, int], bool],
    thr: float,
) -> np.ndarray:
    adjacency_new = np.zeros_like(adjacency)
    for arc in cpdag.arcs:
        i, j = arc
        edge_weight = adjacency[i, j]
        if keep_if(i, j) and edge_weight >= thr:
            adjacency_new[i, j] = edge_weight
    for i, j in cpdag.edges:
        edge_weight = adjacency[i, j] + adjacency[j, i]
        if keep_if(i, j) and edge_weight >= thr:
            adjacency_new[i, j] = adjacency[i, j]
            adjacency_new[j, i] = adjacency[j, i]
    return adjacency_new

def filter_low_confidence_edges(
    adata: AnnData,
    edge_threshold: float,
    *,
    flowsig_network_key: str = "flowsig_network",
    adjacency_key: str = "adjacency",
    filtered_key: str = "filtered",
) -> None:
    """Keep only bootstrap‑stable edges (no biology applied)."""

    net = adata.uns[flowsig_network_key]["network"]
    adjacency = net[adjacency_key]
    cpdag = gpm.PDAG.from_amat(adjacency)

    # simply weight ≥ threshold, no type constraints
    keep = lambda *_: True
    adjacency_filtered = _filter_edges_by_weight(cpdag, adjacency, keep, edge_threshold)

    net[f"{adjacency_key}_{filtered_key}"] = adjacency_filtered
    adata.uns[flowsig_network_key]["network"] = net

def apply_biological_flow(
    adata: AnnData,
    *,
    flowsig_network_key: str = "flowsig_network",
    adjacency_key: str = "adjacency",
    validated_key: str = "validated",
) -> None:
    """
    Enforce inflow → module → (outflow | module)x ordering on a CPDAG.
    Undirected module–module edges are kept bidirectionally.
    """
    flowsig_network = adata.uns[flowsig_network_key]["network"]
    flow_var_info: pd.DataFrame = adata.uns[flowsig_network_key]["flow_var_info"]
    flow_var_types = flow_var_info["Type"].to_numpy()
    adjacency = flowsig_network[adjacency_key]
    cpdag = gpm.PDAG.from_amat(adjacency)

    def biological(i: int, j: int) -> bool:
        type_i, type_j = flow_var_types[i], flow_var_types[j]
        return (
            (type_i == "inflow" and type_j == "module")
            or (type_i == "module" and type_j in {"outflow", "module"})
        )

    validated = _filter_edges_by_weight(cpdag, adjacency, biological, 0.0)
    flowsig_network[f"{adjacency_key}_{validated_key}"] = validated

def construct_intercellular_flow_network(
    adata: AnnData,
    *,
    flowsig_network_key: str = "flowsig_network",
    adjacency_key: str = "adjacency",
) -> nx.DiGraph:
    """
    Convert a (possibly filtered / validated) CPDAG adjacency into a weighted
    directed NetworkX graph, orienting undirected edges by the same biological
    rules as `apply_biological_flow`.
    """
    net = adata.uns[flowsig_network_key]["network"]
    flow_var_info: pd.DataFrame = adata.uns[flowsig_network_key]["flow_var_info"]
    types = flow_var_info["Type"].to_numpy()
    names = flow_var_info.index.to_numpy()
    adjacency = net[adjacency_key]
    totals = _total_edge_weights(adjacency)
    cpdag = gpm.PDAG.from_amat(adjacency)

    flowsig_network = nx.DiGraph()

    def add_edge(i: int, j: int, w: float) -> None:
        flowsig_network.add_edge(names[i], names[j], weight=w)
        flowsig_network.nodes[names[i]]["type"] = types[i]
        flowsig_network.nodes[names[j]]["type"] = types[j]
    
    # Add directed edges
    for i, j in cpdag.arcs:
        if types[i] == "inflow" and types[j] == "module":
            add_edge(i, j, adjacency[i, j] / totals[(min(i, j), max(i, j))])
        elif types[i] == "module" and types[j] in {"outflow", "module"}:
            add_edge(i, j, adjacency[i, j] / totals[(min(i, j), max(i, j))])

    # Add undirected edges
    for i, j in cpdag.edges:
        pair = (min(i, j), max(i, j))
        edge_weight = min(totals[pair], 1.0)
        if types[i] == "inflow" and types[j] == "module":
            add_edge(i, j, edge_weight)
        elif types[i] == "module" and types[j] == "inflow":
            add_edge(j, i, edge_weight)
        elif types[i] == "module" and types[j] == "outflow":
            add_edge(i, j, edge_weight)
        elif types[i] == "outflow" and types[j] == "module":
            add_edge(j, i, edge_weight)
        elif types[i] == types[j] == "module":
            add_edge(i, j, edge_weight)
            add_edge(j, i, edge_weight)

    return flowsig_network