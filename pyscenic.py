import pandas as pd
import numpy as np
import os, glob
import pickle
import pyscenic

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns

if __name__ == '__main__':

    DATA_FOLDER="/u/home/c/cadenmcq/scenic_output"
    RESOURCES_FOLDER="/u/home/c/cadenmcq/scenic_output"
    DATABASE_FOLDER = "/u/home/c/cadenmcq/scenic_output/gene_based"

    DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "mm9-*.mc9nr.genes_vs_motifs.rankings.feather")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")

    MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_tfs.txt')
    SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "PRN_6wk_8wk_matrix.csv")

    ADJACENCIES_FNAME = os.path.join(DATA_FOLDER, "adjacencies.tsv")
    MODULES_FNAME = os.path.join(DATA_FOLDER, "modules.p")
    MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")
    REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")

    N_SAMPLES = 500
    ex_matrix = pd.read_csv(SC_EXP_FNAME, sep=',').T

    tf_names = load_tf_names(MM_TFS_FNAME)

    db_fnames = glob.glob(DATABASES_GLOB)
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

    adjancencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=True)
    adjacencies.to_csv(ADJACENCIES_FNAME, index=False, sep='\t')
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
    with open(MODULES_FNAME, 'wb') as f:
        pickle.dump(modules, f)
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
    df.to_csv(MOTIFS_FNAME)
    regulons = df2regulons(df)
    with open(REGULONS_FNAME, 'wb') as f:
        pickle.dump(regulons, f)