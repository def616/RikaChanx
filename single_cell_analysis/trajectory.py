import scanpy as sc
import seaborn as sns
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colormaps as cm
import numpy as np
from scipy.sparse import issparse
import pandas as pd
import os
import cellrank as cr
import scvelo as scv
from moscot.problems.time import TemporalProblem
from cellrank.kernels import RealTimeKernel

def main():
    import os
    os.environ["CELLRANK_DISABLE_JIT"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["NUMBA_DISABLE_JIT"] = "1"

    import scanpy as sc
    import cellrank as cr
    import scipy.stats as st
    from cellrank.kernels import RealTimeKernel
    from moscot.problems.time import TemporalProblem

    adata = sc.read(file2)

    adata.obs['day'] = adata.obs['sample'].replace({'JJ003': '13.5', 'JJ004': '12.5', 'JJ005': '11.5'})
    adata.obs["day"] = adata.obs["day"].astype(float).astype("category")
    adata.obs["day_numerical"] = adata.obs["day"].astype(float) 

    tp1 = TemporalProblem(adata)
    tp1 = tp1.prepare(time_key="day")
    tp1 = tp1.solve(epsilon=1e-3, tau_a=0.95, scale_cost="mean")
    tmk1 = RealTimeKernel.from_moscot(tp1)
    tmk1.compute_transition_matrix(self_transitions="all", conn_weight=0.2, threshold="auto")

    # anything from here need to be run inside the wrapper function or some issues w/ 
    # multi-processing will show up
    tmk1.plot_random_walks(
        max_iter=300,
        start_ixs={'day': 11.5},
        basis="X_umap",
        seed=0,
        dpi=150,
        size=30,
        n_jobs=1,
    )
    # viewing using mass probability flow
    ax = tmk1.plot_single_flow(
        cluster_key="leiden_0.6",
        time_key="day",
        cluster="2",
        min_flow=0.1,
        xticks_step_size=1,
        show=False,
        ascending=False,
        clusters=None,
    )

    ax = tmk1.plot_single_flow(
        cluster_key="leiden_0.6",
        time_key="day",
        cluster="8",
        min_flow=0.1,
        xticks_step_size=1,
        show=False,
        ascending=False,
        clusters=None,
    )

    ax = tmk1.plot_single_flow(
        cluster_key="leiden_0.6",
        time_key="day",
        cluster="3",
        min_flow=0.1,
        xticks_step_size=1,
        show=False,
        ascending=False,
        time_points=None,
        clusters=None,
    )

    adata.obsp["transition_matrix"] = tmk1.transition_matrix.copy()
    adata.uns["real_time_kernel"] = tmk1
    adata.uns["transition_params"] = tmk1.params

    return adata, tmk1

if __name__ == "__main__":
    adata, tmk1 = main()
    g = cr.estimators.GPCCA(tmk1)
    g.compute_macrostates()
    g.plot_macrostates()

g = cr.estimators.GPCCA(tmk1)

# compute macrostates (groups of cells that behave similarly over time)
g.fit(cluster_key="leiden_0.6", n_states=10) 

# predict terminal macrostates
g.predict_terminal_states()
# plot terminal states
# each cell is colored according to the terminal state it most likely belongs to; higher color intensity reflects greater confidence in the assignment
g.plot_macrostates(which="terminal", legend_loc="right", s=100)
g.plot_macrostates(which="terminal", legend_loc="right", discrete=False)

# estimating fate probabilities
g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False)
g.plot_fate_probabilities(legend_loc = 'right', same_plot=True)

# uncover driver genes
driver_clusters = ['0']
delta_df = g.compute_lineage_drivers(
    lineages=["0"], cluster_key="leiden_0.6", clusters=driver_clusters)
delta_df.head(10)









































