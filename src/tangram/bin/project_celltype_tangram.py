#!/usr/bin/env python3
#
#
#
import argparse
import os
import warnings
import torch
import tangram as tg
import scanpy as sc
import pandas as pd
import numpy as np


### options
parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    type=argparse.FileType('r'),
    help='Input spatial h5ad file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output spatial h5ad file.'
)

parser.add_argument(
    '--mapping',
    type=argparse.FileType('r'),
    dest='map',
    help='Input mapping h5ad file.'
)

parser.add_argument(
    '-r',
    '--reference',
    dest='ref',
    type=argparse.FileType('r'),
    help='Single cell reference h5ad file containing cell type annotation'
)

parser.add_argument(
    '-a',
    '--annotation-celltype',
    dest='anno',
    type=str,
    default='cell.type',
    help="Annotation for selecting from marker genes (default: '%(default)s')"
)


args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]


### main

# I/O
# Expects h5ad file
try:
    adata_spatial = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

# read scRNAseq file
# Expects h5ad file
try:
    adata_ref = sc.read_h5ad(filename=args.ref.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files for refernce single cell.")

# read mapping file
# Expects h5ad file
try:
    adata_map = sc.read_h5ad(filename=args.map.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files for tangram mapping.")

# get mappings
tg.ut.project_cell_annotations(adata_map, adata_spatial, annotation=args.anno)
df_annotations = adata_spatial.obsm["tangram_ct_pred"]

# asssign best cell type
idx = [np.array(df_annotations[df_annotations.index == cellid]).argmax() for cellid in df_annotations.index ]
best_celltypes = df_annotations.columns[idx]
best_celltype_column_name = "tangram_best_" + args.anno


# add cell type annotations use 'n_tangram_' prefix
prefix = "n_tangram_"
df_annotations = df_annotations.add_prefix(prefix)
assert adata_spatial.obs.index.equals(df_annotations.index)

# add to obs
df_annotations.index.name = None
adata_spatial.obs = pd.concat([adata_spatial.obs, df_annotations], axis=1)
adata_spatial.obs[best_celltype_column_name] = best_celltypes

# write output
adata_spatial.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
