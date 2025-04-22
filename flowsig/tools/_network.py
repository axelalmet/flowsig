from typing import Tuple, Optional, Set, Sequence
from collections.abc import Callable
from collections import defaultdict
from dataclasses import dataclass
import numpy as np
from causaldag import unknown_target_igsp, gsp
from causaldag import partial_correlation_suffstat, partial_correlation_test, MemoizedCI_Tester
from causaldag import gauss_invariance_suffstat, gauss_invariance_test, MemoizedInvarianceTester
import time
from ..preprocessing import FlowSigConfig

from tqdm import tqdm
from joblib import Parallel, delayed
from anndata import AnnData
import warnings
warnings.filterwarnings('ignore')

def _resample_indices(
        X: np.ndarray,
        rng: np.random.Generator,
        *,
        indices_by_blocks: Optional[list[np.ndarray]] = None,
) -> np.ndarray:
    n = X.shape[0]
    if indices_by_blocks is None:                       # Sample individual cells
        return X[rng.integers(0, n, n)]

    # Block bootstrapping for spatial data
    X_rs = X.copy()                                
    for indices in indices_by_blocks:
        X_rs[indices] = X[rng.choice(indices, size=indices.size)]
    return X_rs

def _indices_by_block(labels: np.ndarray) -> list[np.ndarray]:
    # inverse gives the block‑id of each row in a single pass
    _, inverse = np.unique(labels, return_inverse=True)
    n_blocks = inverse.max() + 1
    return [np.where(inverse == b)[0] for b in range(n_blocks)]

def _drop_zero_sd_cols(
        matrices: Sequence[np.ndarray]
) -> Tuple[np.ndarray, list[np.ndarray]]:

    # boolean masks of shape (n_features,) for each matrix
    non_zero_masks = [(m.std(0) != 0) for m in matrices]
    keep_mask = np.logical_and.reduce(non_zero_masks)
    keep_idx = np.where(keep_mask)[0]
    filtered = [m[:, keep_idx] for m in matrices]
    return keep_idx, filtered

def _bootstrap_network(
        X_ctrl: np.ndarray,
        X_pert_list: Optional[list[np.ndarray]] = None,
        *,
        indices_by_blocks_ctrl: Optional[np.ndarray] = None,
        indices_by_blocks_pert: Optional[list[np.ndarray]] = None,
        rng: np.random.Generator,
        learner:  Callable[..., Tuple[np.ndarray, ...]],
) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, list[Set[int]]]:
    """
    # Bootstrap samples from data matrix X, either by individual cells or by blocks (for spatial data).
    # Also returns which features have zero standard deviation and should be dropped for instance.
    """
    X_ctrl_rs = _resample_indices(X_ctrl, rng, indices_by_blocks=indices_by_blocks_ctrl)

    if X_pert_list is not None:
        X_pert_rs_list = [
            _resample_indices(X_pert, rng, indices_by_blocks=indices_by_blocks_pert)
            for X_pert, indices_by_blocks in zip(X_pert_list, indices_by_blocks_pert or [None]*len(X_pert_list))
        ]
        keep_idx, mats = _drop_zero_sd_cols([X_ctrl_rs, *X_pert_rs_list])
        X_ctrl_rs, *X_pert_rs_list = mats

        A, pert_targets = learner(X_ctrl_rs, X_pert_rs_list)
        return keep_idx, A, pert_targets
    else:
        # GSP / no perturbed conditions
        keep_idx, (X_ctrl_rs,) = _drop_zero_sd_cols([X_ctrl_rs])
        A = learner(X_ctrl_rs)
        return keep_idx, A

def _learn_gsp(X: np.ndarray, 
               alpha: float = 1e-3
) -> np.ndarray:                      
    suff = partial_correlation_suffstat(X, invert=True)
    ci   = MemoizedCI_Tester(partial_correlation_test, suff, alpha)
    return gsp(set(range(X.shape[1])), ci, nruns=20).cpdag().to_amat()[0]

def _learn_utigsp(X_ctrl: np.ndarray,
                  X_pert_list: list[np.ndarray],
                  alpha: float = 1e-3,
                  alpha_inv: float = 1e-3) -> Tuple[np.ndarray, list[Set[int]]]:
    """
    Learner function for the UT-IGSP algorithm.
    """
    suff_ctrl = partial_correlation_suffstat(X_ctrl, invert=True)
    ci   = MemoizedCI_Tester(partial_correlation_test, suff_ctrl, alpha=alpha)
    
    # Invariance testing
    suff_inv = gauss_invariance_suffstat(X_ctrl, X_pert_list)
    invariance_tester = MemoizedInvarianceTester(gauss_invariance_test, suff_inv, alpha=alpha_inv)

    # Assume unknown interventions for UT-IGSP
    setting_list = [dict(known_interventions=[]) for _ in X_pert_list]

    # Run the UT-IGSP algorithm
    est_dag, est_targets_list = unknown_target_igsp(setting_list, 
                                                    set(range(X_ctrl.shape[1])),
                                                    ci,
                                                    invariance_tester,
                                                    nruns=20)
    
    est_icpdag = est_dag.interventional_cpdag(est_targets_list, cpdag=est_dag.cpdag())

    return est_icpdag.to_amat()[0], est_targets_list

def run_gsp(samples: np.ndarray,
            use_spatial: bool = False,
            indices_by_blocks: Optional[list[np.ndarray]] = None,
            alpha: float = 1e-3,
            seed: int = 0) -> dict:
    
    if use_spatial and indices_by_blocks is None:
        raise ValueError("Block labels must be provided for spatial data.")
    
    # Reseed the random number generator
    rng = np.random.default_rng(seed)
    
    keep, A = _bootstrap_network(X_ctrl=samples,
                                rng=rng,
                                learner=lambda X: _learn_gsp(X, alpha),
                                indices_by_blocks_ctrl=indices_by_blocks if use_spatial else None)
    return {"flow_var_indices": keep, "adjacency_cpdag": A}

def run_utigsp(samples_ctrl: np.ndarray,
               samples_pert_list: list[np.ndarray],
                use_spatial: bool = False,
                indices_by_blocks_ctrl: Optional[list[np.ndarray]] = None,
                indices_by_blocks_pert: Optional[list[np.ndarray]] = None,
                alpha: float=1e-3,
                alpha_inv: float = 1e-3,
                seed: int = 0) -> dict:
    
    if use_spatial and (indices_by_blocks_ctrl is None or indices_by_blocks_pert is None):
        raise ValueError("Block labels must be provided for spatial data (separate for control and perturbed).")
    
    # Reseed the random number generator
    rng = np.random.default_rng(seed)

    keep, A, pert_targets = _bootstrap_network(X_ctrl=samples_ctrl,
                                               X_pert_list=samples_pert_list,
                                               rng=rng,
                                               learner=lambda Xc, Xp: _learn_utigsp(Xc, Xp, alpha, alpha_inv),
                                               indices_by_blocks_ctrl=indices_by_blocks_ctrl if use_spatial else None,
                                               indices_by_blocks_pert=indices_by_blocks_pert if use_spatial else None)
    
    return {'flow_var_indices':keep, 'adjacency_cpdag':A,'perturbed_targets_indices':pert_targets}

# Class to help with bootstrapping
@dataclass(frozen=True, slots=True)
class BootstrapPlan:
    job_id:      int
    seed:        int
    ctrl_X:      np.ndarray
    pert_Xs:     Optional[list[np.ndarray]]
    indices_by_blocks_ctrl: Optional[list[np.ndarray]]
    indices_by_blocks_pert: Optional[list[np.ndarray]]
    learner:     Callable[..., tuple[np.ndarray, ...]]
    alpha_ci:    float 
    alpha_inv:   float

    def run(self) -> dict:
        rng = np.random.default_rng(self.seed)
        if self.pert_Xs is None:
            keep, A = _bootstrap_network(
                self.ctrl_X,
                rng=rng,
                indices_by_blocks_ctrl=self.indices_by_blocks_ctrl,
                learner=lambda X: self.learner(X, self.alpha_ci),
            )
            return {"keep": keep, "A": A}
        else:
            keep, A, pert = _bootstrap_network(
                self.ctrl_X,
                self.pert_Xs,
                rng=rng,
                indices_by_blocks_ctrl=self.indices_by_blocks_ctrl,
                indices_by_blocks_pert=self.indices_by_blocks_pert,
                learner=lambda Xc, Xp: self.learner(
                    Xc, Xp, self.alpha_ci, self.alpha_inv
                ),
            )
            return {"keep": keep, "A": A, "targets": pert}
        
def learn_intercellular_flows(adata: AnnData,
                        condition_key: str | None = None,
                        control: str | None = None, 
                        use_spatial: bool = False,
                        block_key: str | None = None,
                        n_jobs: int = 1,
                        n_bootstraps: int = 100,
                        alpha_ci: float = 1e-3,
                        alpha_inv: float = 1e-3,
                        config: FlowSigConfig = FlowSigConfig(),
) -> None:
    """
    Learn the causal signaling network from cell-type-ligand expression constructed
    from scRNA-seq and a base network derived from cell-cell communication inference.

    This method splits the cell-type-ligand expression into control and perturbed
    samples (one sample for each perturbed condition). We then use UT-IGSP [Squires2020]
    and partial correlation testing to learn the causal signaling DAG and the list of 
    perturbed (intervention) targets.
    
    The base network is also used as a list of initial node permutations for DAG learning.
    To overcome the DAG assumption, as cell-cell communication networks are not necessarily
    DAGS, we use bootstrap aggregation to cover as many possible causal edges and the list
    of node permutations is constructed from all possible DAG subsets of the base network.
    Each boostrap sample is generated by sampling with replacement.

    Parameters
    ----------
    adata
        The annotated dataframe (typically from Scanpy) of the single-cell data.
        Must contain constructed flow expression matrices and knowledge of
        possible cellular flow variables.

    condition_key 
        The label in adata.obs which we use to partition the data.

    control
        The category in adata.obs[condition_key] that specifies which cells belong 
        to the control condition, which is known in causal inference as the observational 
        data.

    flowsig_network_key
        The label for which output will be stored in adata.uns

    flow_expr_key
        The label for which the augmente dflow expression expression is stored in adata.obsm

    use_spatial
        Boolean for whether or not we are analysing spatial data, and thus need to use
        block bootstrapping rather than normal bootstrapping, where we resample across all
        cells.

    block_key
        The label that specfies from which observation key we use to construct (hopefully)
        spatially correlated blocks used for block bootstrapping to learn spatially resolved
        cellular flows. These blocks can be simply just dividing the tissue into rougly
        equally spaced tissue regions, or can be based on tissue annotation (e.g. organ, cell type).
    
    n_jobs
        Number of CPU cores that are used during bootstrap aggregation. If n_jobs > 1, jobs
        are submitted in parallel using multiprocessing

    n_boostraps
        Number of bootstrap samples to generate for causal DAG learning.

    alpha_ci
        The significance level used to test for conditional independence
    
    alpha_inv
        The significance level used to test for conditional invariance.

    Returns
    -------
    flow_vars
        The list of cell-type-ligand pairs used during causal structure learning,
        stored in adata.uns[flowsig_network_key]['flow_vars'].

    adjacency
        The weighted adjacency matrix encoding a bagged CPDAG,
        where weights are determined from bootstrap aggregation. Stored in 
        adata.uns[flowsig_network_key]['adjacency']

    perturbed_targets
        The list of inferred perturbed targets, as determined by conditional invariance
        testing and their bootstrapped probability of perturbations. Stored in
        adata.uns[flowsig_network_key]['perturbed_targets']

    References
    ----------

    .. [Squires2020] Squires, C., Wang, Y., & Uhler, C. (2020, August). Permutation-based
     causal structure learning with unknown intervention targets. In Conference on
     Uncertainty in Artificial Intelligence (pp. 1039-1048). PMLR.

    """

    # Extract the control and perturbed samples
    flow_expr_key = config.flowsig_expr_key
    flowsig_network_key = config.flowsig_network_key

    if flow_expr_key not in adata.obsm.keys():
        raise ValueError(f'flow expression key {flow_expr_key} not found in adata.obsm')
    if flowsig_network_key not in adata.uns.keys():
        raise ValueError(f'flow signature key {flowsig_network_key} not found in adata.uns')
    if use_spatial and block_key is None:   
        raise ValueError("Block key must be provided for spatial data.")
    if block_key and block_key not in adata.obs.keys():
        raise ValueError(f'block key {block_key} not found in adata.obs')
    
    # Get the flow expression matrices
    X_all = adata.obsm[flow_expr_key]
    blocks_all = adata.obs[block_key].values if (use_spatial and block_key) else None
    flow_vars = adata.uns[flowsig_network_key]["flow_var_info"].index.to_list()

    # Use GSP if no perturbation specified
    if condition_key is None:
        ctrl_X, pert_Xs = X_all, None
        ctrl_blocks, pert_blocks = blocks_all, None
        learner = _learn_gsp 
    # Use UT-IGSP for ctrl vs. perturbed
    else:
        if control is None:
            raise ValueError("control must be specified when condition_key is given")

        mask_ctrl = adata.obs[condition_key] == control
        ctrl_X = X_all[mask_ctrl]
        ctrl_blocks= blocks_all[mask_ctrl] if blocks_all is not None else None

        pert_keys = [k for k in adata.obs[condition_key].unique() if k != control]
        pert_Xs = [X_all[adata.obs[condition_key] == k] for k in pert_keys]
        pert_blocks = (
            [blocks_all[adata.obs[condition_key] == k] for k in pert_keys]
            if blocks_all is not None else None
        )
        learner = _learn_utigsp

    indices_by_blocks_ctrl = _indices_by_block(ctrl_blocks) if use_spatial else None
    indices_by_blocks_pert = [_indices_by_block(b) for b in pert_blocks] if (use_spatial and pert_blocks) else None

    root_rng = np.random.default_rng(0)
    seeds = root_rng.integers(0, 2**32 - 1, size=n_bootstraps)

    plans = [
        BootstrapPlan(
            job_id=i,
            seed=int(seeds[i]),
            ctrl_X=ctrl_X,
            pert_Xs=pert_Xs,
            indices_by_blocks_ctrl=indices_by_blocks_ctrl,
            indices_by_blocks_pert=indices_by_blocks_pert,
            learner=learner,
            alpha_ci=alpha_ci,
            alpha_inv=alpha_inv,
        )
        for i in range(n_bootstraps)
    ]

    print(f"Starting {n_bootstraps} bootstraps on {n_jobs} cores …")
    t0 = time.perf_counter()
    with Parallel(n_jobs=n_jobs, prefer="processes") as pool, tqdm(total=n_bootstraps) as bar:
        results = []
        for res in pool(delayed(BootstrapPlan.run)(p) for p in plans):
            results.append(res)
            bar.update()
    elapsed = time.perf_counter() - t0
    print(f"Finished in {elapsed:,.1f} s")

    # Bootstrap aggregation
    bagged_A = np.zeros((len(flow_vars), len(flow_vars)), dtype=float)
    if pert_Xs is not None:
        bagged_targets = defaultdict(lambda: np.zeros(len(flow_vars), dtype=float))

    for res in results:
        keep = res["keep"]
        bagged_A[np.ix_(keep, keep)] += res["A"]
        if "targets" in res:
            for i, tset in enumerate(res["targets"]):
                bagged_targets[i][list(tset)] += 1

    bagged_A /= n_bootstraps
    if pert_Xs is not None:
        bagged_targets = {k: v / n_bootstraps for k, v in bagged_targets.items()}

    # Store output
    adata.uns[flowsig_network_key]["network"] = (
        {
            "flow_vars": flow_vars,
            "adjacency": bagged_A,
            "perturbed_targets": list(bagged_targets.values()),
        }
        if pert_Xs is not None
        else {"flow_vars": flow_vars, "adjacency": bagged_A}
    )