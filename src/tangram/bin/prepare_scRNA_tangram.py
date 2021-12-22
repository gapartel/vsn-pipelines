#!/usr/bin/env python3
#
#
#
import argparse
import os
import warnings
import scanpy as sc
import pandas as pd
import numpy as np

### options
parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input scRNAseq h5ad file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output spatial h5ad file.'
)

parser.add_argument(
    '-m',
    '--method-gene-selection',
    dest='method',
    choices=['marker_genes', 'list', 'all'],
    default='marker_genes',
    help="Method for selecting genes for cell type mapping (default: '%(default)s')"
)

parser.add_argument(
    '-a',
    '--annotation-celltype',
    dest='anno',
    type=str,
    default='cell.type',
    help="Annotation for selecting from marker genes (default: '%(default)s')"
)

parser.add_argument(
    "--normalize",
    dest='do_normalize',
    choices=['true', 'false'],
    default='true',
    help="Normalize data for computing gene rank groups."
)

parser.add_argument(
    '--rank-gene-method',
    dest='method_rank_genes',
    choices=['logreg', 't-test', 'wilcoxon', 't-test_overestim_var'],
    default='wilcoxon',
    help="Method for computing 'rank_genes_groups' in reference data if not present (default: '%(default)s')"
)


args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]


### main

# I/O
# Expects h5ad file
try:
    adata_ref = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

# remove cells and genes with 0 counts everywhere in reference
sc.pp.filter_cells(adata_ref, min_genes=1)
sc.pp.filter_genes(adata_ref, min_cells=1)

# copy expression data
rawX = adata_ref.X.copy()

    
# get marker genes list
gene_list = []

if args.method == 'marker_genes':

    if args.anno not in adata_ref.obs.keys():
        raise Exception("VSN ERROR: Annotation '{}' not found in reference data set.".format(args.anno))
    else:
        if not (hasattr(adata_ref, 'uns') and 'rank_genes_groups' in adata_ref.uns.keys()):
            # compute rank genes groups

            # normalize
            if args.do_normalize == 'true':
                sc.pp.normalize_total(adata_ref)
                
            # get log for ranking genes
            sc.pp.log1p(adata_ref)

            print("Computing 'rank_genes_groups' ...")
            sc.tl.rank_genes_groups(adata_ref, args.anno, method=args.method_rank_genes)
            print("Done.")

            # copy back raw data
            adata_ref.X = rawX.copy()

# to prevent pot. encoding issue substitute tuple entries with strings
for obskey in adata_ref.obs.keys():
    if isinstance(adata_ref.obs[obskey][0], tuple):
        adata_ref.obs[obskey] = [ str(tupl) for tupl in adata_ref.obs[obskey] ]
for varkey in adata_ref.var.keys():
    if isinstance(adata_ref.var[varkey][0], tuple):
        adata_ref.var[varkey] = [ str(tupl) for tupl in adata_ref.var[varkey] ]

# write output
adata_ref.write("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
