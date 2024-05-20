# FlowSig
Python package to infer directed intercellular flows described by ligand-receptor interactions driving tissue-scale gene expression patterning. 

FlowSig requires:

1. Single-cell RNA-sequencing (scRNA-seq) data with cell state annotations that compare a baseline control to one or more perturbed conditions, e.g. disease severity, OR spatial transcriptomics (ST) data.
2. Cell-cell communication (CCC) inferred for each condition of interest. For non-spatial data, we require input from [CellChat](https://github.com/sqjin/CellChat). For ST data, we require input from [COMMOT](https://github.com/zcang/COMMOT).

The code used to generate all of the results in Almet et al., "Inferring pattern-driving intercellular flows from single-cell and spatial transcriptomics", can be found in another GitHub repository [here](https://github.com/axelalmet/FlowSigAnalysis_2023).

<details>
  <summary>Installation</summary>
  
## Installation instructions
  
The easiest way to currently install FlowSig is to generate a Python virtual environment and clone the repository, so that you can install all of the relevant dependencies, particularly those needed by [pyliger](https://github.com/welch-lab/pyliger) and [NSF](https://github.com/willtownes/spatial-factorization-py). We are working on making flowsig pip installable ASAP!

To generate a virtual environment, run the command. N.B. make sure you're using Python 3.8, for some reason pyliger does not like Python 3.9. At least, this is the case on an M1/M2 MacBook.

```
# Create the virtual environment
python3 -m venv flowsigenv

# Activate the virtual environment
source flowsigenv/bin/activate

# Clone the repository
git clone https://github.com/axelalmet/flowsig.git
cd ./flowsig/

# Install using the setup.cfg
pip3 install .
```


You will need to install NSF to infer spatial GEMS separately (see [here](https://github.com/willtownes/spatial-factorization-py))
```
pip3 install git+https://github.com/willtownes/spatial-factorization-py.git#egg=spatial-factorization

```
</details>
<details>
  <summary>Non-spatial scRNA-seq example </summary>

## Application to non-spatial scRNA-seq of stimulated pancreatic islets

Here, we show how to apply FlowSig to an scRNA-seq dataset of wildtype
and stimulated human pancreatic islets, as originally studied in [Burkhardt et al. (2021)](https://www.nature.com/articles/s41587-020-00803-5).
The processed data and cell-cell communication inference, which was obtained using CellChat,
can be downloaded from the following Zenodo [repository](https://zenodo.org/doi/10.5281/zenodo.10850397).

You can also look at the code in a Jupyter notebook found [here](https://github.com/axelalmet/flowsig/blob/main/flowsig/tutorials/pancreatic_islets_scrnaseq_example.ipynb).

### Import packages
```
import flowsig as fs
import scanpy as sc
import pandas as pd
```

### Load the data and cell-cell communication inference

Data is specified in the form of a [Scanpy](https://scanpy.readthedocs.io/en/stable/) object, which is really just an annotated dataframe, i.e. [AnnData](https://anndata.readthedocs.io/en/latest/) object. All subsequent output generated from FlowSig is stored in the Scanpy object.

```
data_directory = '../data/'

# Load the scanpy object
adata = sc.read(data_directory + 'burkhardt21_merged.h5ad')
condition_key = 'Condition'

# Load the cell-cell communication inference
cellchat_Ctrl = pd.read_csv('../communication_inference/output/burkhardt21_leiden_communications_Ctrl.csv')
cellchat_IFNg = pd.read_csv('../communication_inference/output/burkhardt21_leiden_communications_IFNg.csv')

cellchat_output_key = 'cellchat_output'
# Make sure your keys for the cellchat output dictionary match the relevant condition labels
adata.uns[cellchat_output_key] = {'Ctrl': cellchat_Ctrl,
                                  'IFNg': cellchat_IFNg}
```

### Construct GEMs
We now construct gene expression modules (GEMs) from the unnormalised count data. For non-spatial scRNA-seq where we have multiple conditions, we use the iNMF algorithm by [pyliger](https://github.com/welch-lab/pyliger).

```
fs.pp.construct_gems_using_pyliger(adata,
                                n_gems = 10,
                                layer_key = 'counts',
                                condition_key = condition_key)
```

### Construct the flow expression matrices
We construct augmented flow expression matrices for each condition that measure three types of variables:
1. Intercellular signal inflow, i.e., how much of a signal did a cell _receive_. For non-spatial scRNA-seq, signal inflow is defined as receptor gene expression weighted by the average expression of immediate downstream transcription factors that indicate signal activation.
2. GEMs, which encapsulate intracellular information processing. We define these as cellwise membership to the GEM.
3. Intercellular signal outflow, i.e., how much of a signal did a cell _send_. These are simply ligand gene expression.


The kay assumption of flowsig is that all intercellular information flows are directed from signal inflows to GEMs, from one GEM to another GEM, and from GEMs to signal outflows.

For non-spatial scRNA-seq, we need to specify the model organism so that FlowSig knows which receptor-transcription factor targets list to look at.
```
fs.pp.construct_flows_from_cellchat(adata,
                                cellchat_output_key,
                                gem_expr_key = 'X_gem',
                                scale_gem_expr = True,
                                model_organism = 'human',
                                flowsig_network_key = 'flowsig_network',
                                flowsig_expr_key = 'X_flow')
```


To reduce the number of variables over which we have to infer intercellular flows—and thus computation time—and to prioritise 'informative variables', we only retain inflow and outflow signals that are sufficiently _differentially flowing_ between the control and perturbed conditions. We determine differentially flowing signals using a Wilcoxon rank-sum test and retain variables only if they are below a specified adjusted p-value threshold (q-value) and above a specified log-fold-change threshold.

```
fs.pp.determine_informative_variables(adata,  
                                    flowsig_expr_key = 'X_flow',
                                    flowsig_network_key = 'flowsig_network',
                                    spatial = False,
                                    condition_key = condition_key,
                                    control_key =  'Ctrl',
                                    qval_threshold = 0.05,
                                    logfc_threshold = 0.5)
```

### Learn intercellular flows

We are now in a position to learn the intercellular flows. To increase reliability of objects, we bootstrap aggregate results over a number of realisations. For non-spatial data, we have to specify the condition label and the control condition.

This step uses [UT-IGSP](https://uhlerlab.github.io/causaldag/utigsp.html) to learn what is called a completed partially directed acyclic graph (CPDAG), which encodes directed arcs and undirected edges that describe the Markov Equivalence Class of statistical dependence relations that were learned directly from the data using conditional independence testing (how do variables depend on one another) and conditional invariance testing (which variables changed significantly between conditions). For both tests, we use a parametric partial-correlation-based method. The main reason we used these tests were because they take the least time to run compared to nonparametric kernel-based tests. Any test like the Hilbert-Schmidt Independence Criterion takes way too long for even 10-20 variables. The big caveat is that partial correlation assumes the data is described by a linear Gaussian model, which obviously isn't true for scRNA-seq. It's a long-term goal to add different types of nonparametric conditional independence/invariance tests that can be run in a reasonable amount of time. 

```
fs.tl.learn_intercellular_flows(adata,
                        condition_key = condition_key,
                        control_key = 'Ctrl', 
                        flowsig_key = 'flowsig_network',
                        flow_expr_key = 'X_flow',
                        use_spatial = False,
                        n_jobs = 4,
                        n_bootstraps = 500)
```
### Partially validate intercellular flow network

Finally, we will remove any "false positive" edges. Noting that the CPDAG contains directed arcs and undirected arcs we do two things. 

First, we remove directed arcs that are not oriented from signal inflow to GEM, GEM to GEM, or from GEM to signal outflow and for undirected edges, we reorient them so that they obey the previous directionalities.

```
fs.tl.apply_biological_flow(adata,
                            flowsig_network_key = 'flowsig_network',
                            adjacency_key = 'adjacency',
                            validated_key = 'validated')
```

Second, we will remove directed arcs whose bootstrapped frequencies are below a specified edge threshold as well as undirected edges whose total bootstrapped frequencies are below the same threshold.

```
edge_threshold = 0.7
fs.tl.filter_low_confidence_edges(adata,
                                edge_threshold = edge_threshold,
                                flowsig_network_key = 'flowsig_network',
                                adjacency_key = 'adjacency_validated',
                                filtered_key = 'filtered')
```

Every time we apply these steps, we generate a new adjacency matrix that describes the intercellular flow network edges. The original output from `learn_intercellular_flows` is stored in `adata.uns['flowsig_network]['network']['adjacency']`. After validation using `apply_biological_flow`, we add the `_validated` key to the adjacency key, so that the new adjacency is stored in `adata.uns['flowsig_network]['network']['adjacency_validated']`. After filtering low-confidence edges, the adjacency is stored under `adata.uns['flowsig_network]['network']['adjacency_validated_filtered']`. As a note, you could change the order of these two post-processing steps, which would mean that the final adjacency would be under the key `adjacency_filtered_validated`. The results should be similar but it's worth remembering this.

We also note that if you want to explore the network directly, we included a function to generate the directed [NetworkX](https://networkx.org/documentation/stable/index.html) `DiGraph` object. You will need to generate this for any of the plotting functions.

```
flow_network = fs.tl.construct_intercellular_flow_network(adata,
                                                        flowsig_network_key = 'flowsig_network',
                                                        adjacency_key = 'adjacency_validated_filtered')
```
</details>

<details>
  <summary>Spatial example</summary>
  
  ## Application to spatial Stereo-seq of E9.5 mouse embryo 
Here, we show how to apply FlowSig to a spatial Stereo-seq dataset of an E9.5 mouse embryo, as originally studied in [Chen et al. (2022)](https://doi.org/10.1016/j.cell.2022.04.003).
The processed data and cell-cell communication inference, which was obtained using [COMMOT](https://commot.readthedocs.io/en/latest/tutorials.html),
can be downloaded from the following Zenodo  [repository](https://zenodo.org/doi/10.5281/zenodo.10850397).

You can also look at the code in a Jupyter notebook found [here](https://github.com/axelalmet/flowsig/blob/main/flowsig/tutorials/mouse_embryo_stereoseq_example.ipynb).

### Import packages
```
import flowsig as fs
import scanpy as sc
import pandas as pd
```

### Load the data and cell-cell communication inference

We load the data as an `AnnData` object, which has been subsetted for spatially variable genes only and includes the output from COMMOT already. We note here that COMMOT uses the CellChat database by default and we need to specify where it's been stored.

```
data_directory = '../data/'

# Load the scanpy object
adata = sc.read(data_directory + 'chen22_E9.5_svg.h5ad')
commot_output_key = 'commot-cellchat'
```

### Construct GEMs
We now construct gene expression modules (GEMs) from the unnormalised count data. For ST data, we use [NSF](https://github.com/willtownes/spatial-factorization-py).

```
fs.pp.construct_gems_using_nsf(adata,
                            n_gems = 20,
                            layer_key = 'count',
                            n_inducing_pts = 500,
                            length_scale = 10)
```

### Construct the flow expression matrices
We construct augmented flow expression matrices for each condition that measure three types of variables:
1. Intercellular signal inflow, i.e., how much of a signal did a cell _receive_. For ST data, signal inflow is constructed by summing the received signals for each significant ligand inferred by COMMOT.
2. GEMs, which encapsulate intracellular information processing. We define these as cellwise membership to the GEM.
3. Intercellular signal outflow, i.e., how much of a signal did a cell _send_. These are simply ligand gene expression.

The kay assumption of flowsig is that all intercellular information flows are directed from signal inflows to GEMs, from one GEM to another GEM, and from GEMs to signal outflows.

For spatial data, we use COMMOT output directly to construct signal inflow expression and do not need knowledge about TF databases.
```
fs.pp.construct_flows_from_commot(adata,
                                commot_output_key,
                                gem_expr_key = 'X_gem',
                                scale_gem_expr = True,
                                flowsig_network_key = 'flowsig_network',
                                flowsig_expr_key = 'X_flow')
```

For spatial data, we retain spatially informative variables, which we determine by calculating the Moran's I value for signal inflow and signal outflow variables. In case the spatial graph has not been calculated for this data yet, FlowSig will do so, meaning that we need to specify both the coordinate type, `grid` or `generic`, and in the case of the former, `n_neighs`, which in this case, is 8.

Flow expression variables are defined to be spatially informative if their Moran's I value is above a specified threshold.

```
fs.pp.determine_informative_variables(adata,  
                                    flowsig_expr_key = 'X_flow',
                                    flowsig_network_key = 'flowsig_network',
                                    spatial = True,
                                    moran_threshold = 0.15,
                                    coord_type = 'grid',
                                    n_neighbours = 8,
                                    library_key = None)
```

### Learn intercellular flows
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
                        flowsig_key = 'flowsig_network',
                        flow_expr_key = 'X_flow',
                        use_spatial = True,
                        block_key = 'spatial_kmeans',
                        n_jobs = 4,
                        n_bootstraps = 500)
```
### Partially validate intercellular flow network

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
</details>
