# Application to spatial Stereo-seq of E9.5 mouse embryo 
Here, we show how to apply FlowSig to a spatial Stereo-seq dataset of an E9.5 mouse embryo, as originally studied in [Chen et al. (2022)](https://doi.org/10.1016/j.cell.2022.04.003).
The processed data and cell-cell communication inference, which was obtained using [COMMOT](https://commot.readthedocs.io/en/latest/tutorials.html),
can be downloaded from the following Zenodo  [repository](https://zenodo.org/doi/10.5281/zenodo.10850397).

You can also look at the code in a Jupyter notebook found [here](https://github.com/axelalmet/flowsig/blob/main/flowsig/tutorials/mouse_embryo_stereoseq_example.ipynb).

## Import packages
```
import flowsig as fs
import scanpy as sc
import pandas as pd
```

## Load the data and cell-cell communication inference

We load the data as an `AnnData` object, which has been subsetted for spatially variable genes only and includes the output from COMMOT already. We note here that COMMOT uses the CellChat database by default and we need to specify where it's been stored.

```
data_directory = '../data/'

# Load the scanpy object
adata = sc.read(data_directory + 'chen22_E9.5_svg.h5ad')
commot_output_key = 'commot-cellchat'
```

## Construct GEMs
We now construct gene expression modules (GEMs) from the unnormalised count data. For ST data, we use [NSF](https://github.com/willtownes/spatial-factorization-py).

```
fs.pp.construct_gems_using_nsf(adata,
                            n_gems = 20,
                            layer_key = 'count',
                            n_inducing_pts = 500)
```

## Construct the flow expression matrices
We construct augmented flow expression matrices for each condition that measure three types of variables:
1. Intercellular signal inflow, i.e., how much of a signal did a cell _receive_. For ST data, signal inflow is constructed by summing the received signals for each significant ligand inferred by COMMOT.
2. GEMs, which encapsulate intracellular information processing. We define these as cellwise membership to the GEM.
3. Intercellular signal outflow, i.e., how much of a signal did a cell _send_. These are simply ligand gene expression.

The kay assumption of flowsig is that all intercellular information flows are directed from signal inflows to GEMs, from one GEM to another GEM, and from GEMs to signal outflows.

For spatial data, we use COMMOT output directly to construct signal inflow expression and do not need knowledge about TF databases.
```
fs.pp.construct_flow_expressions(adata,
                      commot_output_key=commot_output_key,
                      spatial = True)
```

For spatial data, we retain spatially informative variables, which we determine by calculating the Moran's I value for signal inflow and signal outflow variables. In case the spatial graph has not been calculated for this data yet, FlowSig will do so, meaning that we need to specify both the coordinate type, `grid` or `generic`, and in the case of the former, `n_neighs`, which in this case, is 8.

Flow expression variables are defined to be spatially informative if their Moran's I value is above a specified threshold.

```
fs.pp.determine_informative_variables(adata,  
                                    spatial = True,
                                    moran_threshold = 0.15,
                                    coord_type = 'grid',
                                    n_neighbours = 8,
                                    library_key = None)
```

## Learn intercellular flows
For spatial data, where there are far fewer "control _vs._ perturbed" studies, we use the [GSP](https://graphical-model-learning.readthedocs.io/en/latest/dags/generated/graphical_model_learning.gsp.html) method, which uses conditional independence testing and a greedy algorithm to learn the CPDAG containing directed arcs and undirected edges.

For spatial data, we cannot bootstrap by resampling across individual cells because we would lose the additional layer of correlation contained in the spatial data. Rather, we divide the tissue up into spatial "blocks" and resample within blocks. This is known as block bootstrapping.

To calculate the blocks, we used scikit-learn's [k-means](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html) clustering method to generate 20 roughly equally sized spatial blocks.

```
from sklearn.cluster import KMeans

kmeans = KMeans(n_clusters=20, random_state=0).fit(adata.obsm['spatial'])
adata.obs['spatial_kmeans'] = pd.Series(kmeans.labels_, dtype='category').values
```
We use these blocks to learn the spatial intercellular flows.

```
fs.tl.learn_intercellular_flows(adata,
                        use_spatial = True,
                        block_key = 'spatial_kmeans',
                        n_jobs = 4,
                        n_bootstraps = 500)
```
## Partially validate intercellular flow network

Finally, we will remove any "false positive" edges. Noting that the CPDAG contains directed arcs and undirected arcs we do two things. 

First, we remove directed arcs that are not oriented from signal inflow to GEM, GEM to GEM, or from GEM to signal outflow and for undirected edges, we reorient them so that they obey the previous directionalities.

```
fs.tl.apply_biological_flow(adata,
                            flowsig_network_key = 'flowsig_network',
                            adjacency_key = 'adjacency',
                            validated_key = 'validated')
```

Second, we will remove directed arcs whose bootstrapped frequencies are below a specified edge threshold as well as undirected edges whose total bootstrapped frequencies are below the same threshold. Because we did not have perturbation data, we specify a more stringent edge threshold.

```
edge_threshold = 0.8
fs.tl.filter_low_confidence_edges(adata,
                                edge_threshold = edge_threshold,
                                flowsig_network_key = 'flowsig_network',
                                adjacency_key = 'adjacency_validated',
                                filtered_key = 'filtered')
```

We can construct the directed [NetworkX](https://networkx.org/documentation/stable/index.html) `DiGraph` object from `adjacency_validated_filtered`.

```
flow_network = fs.tl.construct_intercellular_flow_network(adata,
                                                        flowsig_network_key = 'flowsig_network',
                                                        adjacency_key = 'adjacency_validated_filtered')
```