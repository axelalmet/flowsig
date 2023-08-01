# flowsig
Python package to infer directed intercellular flows described by ligand-receptor interactions driving tissue-scale gene expression patterning. 

FlowSig requires:

1. Single-cell RNA-sequencing (scRNA-seq) data with cell state annotations that compare a baseline control to one or more perturbed conditions, e.g. disease severity, OR spatial transcriptomics (ST) data.
2. Cell-cell communication (CCC) inferred for each condition of interest. For non-spatial data, we recommend [CellChat](https://github.com/sqjin/CellChat), but input from [CellPhoneDB](https://www.cellphonedb.org/) applied to human data is supported. For ST data, we recommend [COMMOT](https://github.com/zcang/COMMOT).

## Examples:
[Non-spatial scRNA-seq example](#non-spatial)
[ST example](#spatial)

<a name="non-spatial"/>
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

### Construct the flow expression matrices.


### Partially validate causal network

Finally, to remove any "false positive" causal relations that are not related
via signalling mechanisms, we benchmark the learned network against the original base network constructed from CCC inference.


```
# Validate causal network against base network
fs.tl.validate_against_base_network(adata_burkhardt,
                            condition_label='Condition',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks')
```

<a name="spatial"/>
## Application to Spatial Stereo-seq of E9.5 mouse embryo

