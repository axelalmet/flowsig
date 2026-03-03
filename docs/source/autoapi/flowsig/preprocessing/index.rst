flowsig.preprocessing
=====================

.. py:module:: flowsig.preprocessing


Functions
---------

.. autoapisummary::

   flowsig.preprocessing.construct_flows_from_cellchat
   flowsig.preprocessing.construct_flows_from_commot
   flowsig.preprocessing.construct_flows_from_cellphonedb
   flowsig.preprocessing.construct_flows_from_liana
   flowsig.preprocessing.construct_flow_expressions
   flowsig.preprocessing.subset_for_flow_type
   flowsig.preprocessing.filter_flow_vars
   flowsig.preprocessing.determine_differentially_flowing_vars
   flowsig.preprocessing.determine_spatially_flowing_vars
   flowsig.preprocessing.determine_informative_variables
   flowsig.preprocessing.construct_gems_using_pyliger
   flowsig.preprocessing.construct_gems_using_nsf
   flowsig.preprocessing.construct_gems_using_nmf
   flowsig.preprocessing.construct_gems_using_cnmf
   flowsig.preprocessing.construct_spatial_blocks


Package Contents
----------------

.. py:function:: construct_flows_from_cellchat(adata: anndata.AnnData, cellchat_output_key: str, model_organism: Literal['human', 'mouse'] = 'human', tfs_to_use: Optional[list[str]] = None, config: FlowSigConfig = FlowSigConfig(), construction: Literal['v1', 'v2'] = 'v1') -> None

.. py:function:: construct_flows_from_commot(adata: anndata.AnnData, commot_output_key: str, config: FlowSigConfig = FlowSigConfig()) -> None

.. py:function:: construct_flows_from_cellphonedb(adata: anndata.AnnData, cellphonedb_output_key: str, cellphonedb_tfs_key: str, model_organism: str = 'human', config: FlowSigConfig = FlowSigConfig())

.. py:function:: construct_flows_from_liana(adata: anndata.AnnData, liana_output_key: str, use_tfs: bool = False, model_organism: Literal['human', 'mouse'] = 'human', config: FlowSigConfig = FlowSigConfig()) -> None

.. py:function:: construct_flow_expressions(adata: anndata.AnnData, *, spatial: bool, method: Literal['cellchat', 'cellphonedb', 'liana'] | None = None, **kwargs: Any) -> None

.. py:function:: subset_for_flow_type(adata: anndata.AnnData, var_type: str = 'all', config: flowsig.preprocessing.FlowSigConfig = FlowSigConfig()) -> anndata.AnnData

.. py:function:: filter_flow_vars(adata: anndata.AnnData, vars_subset: Sequence[str], config: flowsig.preprocessing.FlowSigConfig = FlowSigConfig()) -> None

.. py:function:: determine_differentially_flowing_vars(adata: anndata.AnnData, condition_key: str, control: str, *, config: flowsig.preprocessing.FlowSigConfig = FlowSigConfig(), logfc_threshold: float = None, qval_threshold: float = None, construction: Literal['v1', 'v2'] = 'v1') -> None

.. py:function:: determine_spatially_flowing_vars(adata: anndata.AnnData, *, config: flowsig.preprocessing.FlowSigConfig = FlowSigConfig(), n_perms: int | None = None, n_jobs: int | None = None, library_key: str | None = None, moran_threshold: float = 0.1, coord_type: str | None = None, n_neighs: int | None = None, n_rings: int | None = None, backend: str = 'loky') -> None

.. py:function:: determine_informative_variables(adata: anndata.AnnData, *, config: flowsig.preprocessing.FlowSigConfig = FlowSigConfig(), spatial: bool = False, condition_key: str | None = None, control: str | None = None, **kwargs) -> None

   Wrapper that chooses spatial vs. differential selection.


.. py:function:: construct_gems_using_pyliger(adata: anndata.AnnData, n_gems: int, layer_key: str, condition_key: str)

.. py:function:: construct_gems_using_nsf(adata: anndata.AnnData, n_gems: int, layer_key: str, spatial_key: str = 'spatial', sample_key: Optional[str] = None, n_inducing_pts: int = 100, weight_prior: str = 'Horseshoe', factor_prior: str = 'GP', likelihood: str = 'NegativeBinomial', kernel: str = 'Matern')

.. py:function:: construct_gems_using_nmf(adata: anndata.AnnData, n_gems: int, layer_key: str, random_state: int = 0, max_iter: int = 1000)

.. py:function:: construct_gems_using_cnmf(adata: anndata.AnnData, n_gems: int, usage_norms: numpy.ndarray | pandas.DataFrame, spectra_scores: numpy.ndarray | pandas.DataFrame, spectra_tpm: numpy.ndarray | pandas.DataFrame, cnmf_vars: Optional[list] = None)

.. py:function:: construct_spatial_blocks(adata: anndata.AnnData, n_blocks: int, use_graph: bool = False, graph_adjacency: str = 'spatial_connectivities', resolution: float = None, spatial_block_key: str = 'spatial_block', spatial_key: str = 'spatial')

