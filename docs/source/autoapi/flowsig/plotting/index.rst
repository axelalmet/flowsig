flowsig.plotting
================

.. py:module:: flowsig.plotting


Functions
---------

.. autoapisummary::

   flowsig.plotting.plot_differentially_flowing_signals
   flowsig.plotting.plot_intercellular_flows


Package Contents
----------------

.. py:function:: plot_differentially_flowing_signals(adata: scanpy.AnnData, condition_key: str, pert_key: Optional[Iterable[str]] = None, var_type: str = 'all', flowsig_expr_key: str = 'X_flow', flowsig_network_key: str = 'flowsig_network', qval_threshold: float = 0.05, logfc_threshold: float = 0.5, label_lowqval: bool = False, scatter_size: float = 50, ax: matplotlib.axes.Axes = None)

.. py:function:: plot_intercellular_flows(adata: scanpy.AnnData, flow_network: networkx.DiGraph, inflow_vars: Union[str, Sequence[str]] = None, module_vars: Union[str, Sequence[str]] = None, outflow_vars: Union[str, Sequence[str]] = None, flowsig_network_key: str = 'flowsig_network', align_mode: str = 'horizontal', width_scale: float = 2.0, x_margin_offset: float = 0.3, y_margin_offset: float = 0.0, ax: matplotlib.axes.Axes = None)

