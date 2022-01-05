#!/usr/bin/env python3
#
#
#
import argparse
import os
import warnings
import SpatialDE
import NaiveDE
import scanpy as sc
import pandas as pd
import numpy as np
import statsmodels.stats.multitest as multitest
import scipy

#########################################################
# helper function: test whether spatialDE q-value and fraction of variance (fsv) are satisfying specified trhesholds
#########################################################
def test_spatialde_signif(qval, fsv, thr_qval, min_fsv):
    
    if qval < thr_qval and fsv >= min_fsv:
        return True
    else:
        return False
    

# options
parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input h5ad file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output h5ad file.'
)

parser.add_argument(
    "--min-cells",
    type=int,
    default=10,
    dest='min_cells',
    help='Filter genes by min. number of cells expressing gene.'
)

parser.add_argument(
    "--thr-qval",
    type=float,
    default=0.05,
    dest='thr_qval',
    help='Q-value threshold for spatialDE results to be considered significant.'
)

parser.add_argument(
    '--method-pval-correction',
    dest='method_qval',
    choices=['qval', 'bonferroni', 'holm', 'fdr_bh'],
    default='bonferroni',
    help="Method for multiple testing used for finding significant genes."
)

parser.add_argument(
    "--min-fsv",
    type=float,
    default=0.5,
    dest='min_fsv',
    help='Minimum fraction of spatial variance for spatialDE results to be considered significant.'
)


parser.add_argument(
    "--normalize-naivede",
    choices=['true', 'false'],
    default='true',
    dest='run_naivede',
    help='Normalize counts using NaiveDE and regress out n_counts.'
)



args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

min_cells = args.min_cells
min_fsv = args.min_fsv
thr_qval = args.thr_qval
method_qval = args.method_qval
run_naivede = args.run_naivede

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")


### main
# read anndata
adata = sc.read_h5ad(filename=FILE_PATH_IN.name)

# filter genes
if min_cells > 0:
    sc.pp.filter_genes(adata, min_cells=min_cells)

# get matrix for spatialIDE
if isinstance(adata.X, scipy.sparse.csr.csr_matrix):
    counts = pd.DataFrame(adata.X.todense(), columns=adata.var_names, index=adata.obs_names)
else:
    counts = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs_names)
coord = pd.DataFrame(adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=adata.obs_names)
sample_info = adata.obs


# normalize data
if run_naivede == 'true':
    print("Running NaiveDE ...")
    norm_expr = NaiveDE.stabilize(counts.T).T
    regress_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(n_counts)').T
    counts = regress_expr
    print("Done.")
    
# run spatialDE
print("Running SpatialDE ...")
results = SpatialDE.run(coord, counts)
print("Done.")

# run additional corrections for multiple testing
for test_method in ['bonferroni', 'holm', 'fdr_bh']:
    results_multitest = multitest.multipletests(results.pval, method=test_method, alpha=thr_qval)
    results[test_method] = results_multitest[1]
    
# check significance
dict_signif = {'is_signif': [test_spatialde_signif(qval, fsv, thr_qval, min_fsv) for qval, fsv in zip(results[method_qval], results['FSV'])]}
results['is_signif'] = dict_signif['is_signif']

# add results to uns
adata.uns['spatialDE'] = results
# update X with (normalized) counts
adata.X = scipy.sparse.csr_matrix(counts)

# write output
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
