import os
import pandas as pd
import scanpy as sc
import seaborn as sns
import pandas as pd
import numpy as np
import glob,os
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
path = r'~'
file = glob.glob(os.path.join(path, "*.csv"))#
print(file)
dl = []
for f in file:
    if os.path.splitext(f)[1] == '.csv':
            fileName = os.path.splitext(f)[0][62:73] #
            adata=sc.read_csv(f).T#
            adata.obs['sample']=fileName
            dl.append(adata)
adata = dl[0].concatenate(dl[1],dl[2],dl[3],dl[4],dl[5],dl[6],dl[7],dl[8])
adata2 = sc.AnnData(X=adata.X, var=adata.var, obs = adata.obs)
sc.pp.normalize_per_cell(adata2, counts_per_cell_after=1e4)
sc.pp.log1p(adata2)
sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata2)
print("Highly variable genes: %d"%sum(adata2.var.highly_variable))
var_genes_all = adata2.var.highly_variable
sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'batch')
print("Highly variable genes intersection: %d"%sum(adata2.var.highly_variable_intersection))
print("Number of batches where gene is variable:")
print(adata2.var.highly_variable_nbatches.value_counts())
var_genes_batch = adata2.var.highly_variable_nbatches > 0
var_select = adata2.var.highly_variable_nbatches > 1
var_genes = var_select.index[var_select]
len(var_genes)
batches = ['0','1','2','3','4','5','6','7','8']
alldata = {}
for batch in batches:
    alldata[batch] = adata2[adata2.obs['batch'] == batch,]

alldata  
import scanorama
alldata2 = dict()
for ds in alldata.keys():
    print(ds)
    alldata2[ds] = alldata[ds][:,var_genes]

#convert to list of AnnData objects
adatas = list(alldata2.values())

# run scanorama.integrate
scanorama  = scanorama.integrate_scanpy(adatas, dimred = 50,)
all_s = np.concatenate(scanorama)
print(all_s.shape)
# add to the AnnData object
adata2.obsm["SC"] = all_s
sc.pp.neighbors(adata2, n_pcs =50, use_rep = "SC")
sc.tl.umap(adata2)
fig, axs = plt.subplots(1, figsize=(11,7),constrained_layout=True)
sc.pl.umap(adata, color="sample", title="scanorama", show=False,size=20, legend_loc='right margin',ax=axs,sort_order=True)