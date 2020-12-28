import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar
#improtant add
from distributed import Client, LocalCluster

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.export import export2loom

import seaborn as sns


#path
if __name__ == '__main__':
    RESOURCES_FOLDER = "~"
    DATA_FOLDER = "~"
    SCHEDULER = "123.122.8.24:8786"
    DATABASES_GLOB = os.path.join(RESOURCES_FOLDER,  "mm9-*.mc9nr.feather")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
    MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_tfs.txt')
    SC_EXP_FNAME = os.path.join(DATA_FOLDER, "counts.txt")
    REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
    MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")

    # preliminary work
    ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
    #ex_matrix3 = ex_matrix.iloc[:, :].values
    tf_names = load_tf_names(MM_TFS_FNAME)

    db_fnames = glob.glob(DATABASES_GLOB)

    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]


    #assert ex_matrix2.shape[1] == len(gene_names)
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

    # Phase I: Inference of co-expression modules
    adjacencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=True)
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))


    # Phase II: Prune modules for targets with cis regulatory footprints
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)

    # Save the enriched motifs and the discovered regulons to disk.
    df.to_csv(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)

    #Phase III: Cellular regulon enrichment matrix (aka AUCell)
    auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
    AUC_FNAME = os.path.join(DATA_FOLDER, "macrophage.csv")
    AUC_FNAME1 = os.path.join(DATA_FOLDER, "macrophage1.csv")
    auc_mtx.to_csv(AUC_FNAME)
    regulons = [r.rename(r.name.replace('(+)', ' (' + str(len(r)) + 'g)')) for r in regulons]
    auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
    auc_mtx.to_csv(AUC_FNAME1)