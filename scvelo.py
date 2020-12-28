import os
os.chdir('~')
import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
from IPython.display import Image
import scanpy as sc
import scvelo as scv
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo')  # for beautified visualization
adata = scv.read('sample1.loom', cache=True)
adata.var_names_make_unique()
scv.pp.filter_and_normalize(adata, min_shared_counts=10, n_top_genes=3000)
scv.pp.moments(adata,n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
sc.tl.umap(adata)
scv.tl.recover_dynamics(adata)

scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata,basis='umap', color='cluster',size=180,alpha=0.8,legend_loc='none')
scv.tl.terminal_states(adata)
scv.pl.scatter(adata, color=[ 'root_cells', 'end_points'],size=80)
