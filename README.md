# flowsig
Python package to infer directed intercellular flows described by ligand-receptor interactions driving tissue-scale gene expression patterning. 

FlowSig requires:

1. Single-cell RNA-sequencing (scRNA-seq) data with cell state annotations that compare a baseline control to one or more perturbed conditions, e.g. disease severity, OR spatial transcriptomics (ST) data.
2. Cell-cell communication (CCC) inferred for each condition of interest. For non-spatial data, we recommend [CellChat](https://github.com/sqjin/CellChat), but input from [CellPhoneDB](https://www.cellphonedb.org/) applied to human data is supported. For ST data, we recommend [COMMOT](https://github.com/zcang/COMMOT).

Examples:
[Non-spatial scRNA-seq example](## Application to non-spatial scRNA-seq of stimulated pancreatic islets)
[ST example](## Application to spatial Stereo-seq of E9.5 mouse embryo)


## Application to non-spatial scRNA-seq of stimulated pancreatic islets
Here, we show how to apply FlowSig to an scRNA-seq dataset of wildtype
and stimulated human pancreatic islets, as originally studied in [Burkhardt et al. (2021)](https://www.nature.com/articles/s41587-020-00803-5).
The processed data and cell-cell communication inference, which was obtained using CellChat,
can be downloaded from the following Zenodo [repository](https://zenodo.org/record/6791333#.Ys80tuxKirc).

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
adata_burkhardt = sc.read(data_directory + 'burkhardt21_merged.h5ad')

# Load the cell-cell communication inference
cellchat_burkhardt_ctrl = pd.read_csv(data_directory + 'CellCellCommunication/burkhardt21_communications_Ctrl.csv')
cellchat_burkhardt_ifng = pd.read_csv(data_directory + 'CellCellCommunication/burkhardt21_communications_IFNg.csv')

```
This next step is entirely optional, but we removed weaker interactions to reduce
noise during causal learning.
```
pval_threshold = 0.05
quantile = 0.25
cellchat_burkhardt_ctrl = cellchat_burkhardt_ctrl[(cellchat_burkhardt_ctrl['pval'] < pval_threshold)&(cellchat_burkhardt_ctrl['prob'] > cellchat_burkhardt_ctrl['prob'].quantile(quantile))]
cellchat_burkhardt_ifng = cellchat_burkhardt_ifng[(cellchat_burkhardt_ifng['pval'] < pval_threshold)&(cellchat_burkhardt_ifng['prob'] > cellchat_burkhardt_ifng['prob'].quantile(quantile))]
```
We define CCC output as a dictionary, where the keys are the conditions to which the
inferred lists of interactions correspond.

```
cellchat_output_burkhardt = {'Ctrl':cellchat_burkhardt_ctrl, 'IFNg':cellchat_burkhardt_ifng}
```

### Construct the base network
We now construct the base network of candidate causal signalling interactions. We construct the network by tracing through both lists of interactions for communicating triplets of cell types, A -> B -> C such that A -> B via one ligand-receptor pair, L1 ~ R1, and B -> C via another ligand-receptor pair, L2 ~ R2. Note that it is possible that A = C and that L1 = L2, R1 = R2.

We also point out that you need to specify which cell-cell communication inference method you used, i.e. cellchat, cellphonedb, or squidpy, and you also need to know the observation column that specifies the condition in your Scanpy object.

```
# Construct the base network now
cfs.pp.construct_base_networks(adata_burkhardt,
                                cellchat_output_burkhardt,
                                condition_label='Condition',
                                method='cellchat',
                                node_sep='_',
                                base_network_label='base_networks')
```
### Construct cell-type-ligand expressions
Rather than perform inference on all ~20K genes, we only consider
all possible cell-type-ligand pairs that were implicated by the base network.
Therefore, we need to transform the scRNA-seq data from all genes to
cell-type-specific ligand expression. This is also what allows us to infer
causal signalling at the cell type resolution.

```
# Construct the cell-type ligand expression matrix
cfs.pp.construct_celltype_ligand_expressions(adata_burkhardt,
                                            celltype_label='Type',
                                            node_sep='_',
                                            expressions_label='X_celltype_ligand',
                                            base_network_label='base_networks')
```
### Learn causal network
We now learn the causal network. For this step, the control condition needs to be prescribed. The results are bootstrap aggregated to improve confidence in the inferred edges, where `n_bootstraps` sets the number of bootstrap samples to genreate.

The computational time to run this step scales up 
significantly as the number of variables increases—as is the case for all causal structure learning methods—so we recommend running this step across multiple cores by increasing the `n_cores` option.

```
# Learn the causal signaling networks now
cfs.tl.learn_causal_network(adata_burkhardt,
                            condition_label='Condition',
                            control_label='Ctrl',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks', 
                            n_cores=4,
                            n_bootstraps=100,
                            n_shuffles=1000,
                            alpha_ci=0.001,
                            alpha_inv=0.001)
```

### Partially validate causal network

Finally, to remove any "false positive" causal relations that are not related
via signalling mechanisms, we benchmark the learned network against the original base network constructed from CCC inference.


```
# Validate causal network against base network
cfs.tl.validate_against_base_network(adata_burkhardt,
                            condition_label='Condition',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks')
```

## Application to Spatial 

