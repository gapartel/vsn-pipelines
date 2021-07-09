#!/usr/bin/env python3
#
#
#
import argparse
import os
import warnings
import SpatialDE
import scanpy as sc
import pandas as pd

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

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")


### main
# read anndata
adata = sc.read_h5ad(filename=FILE_PATH_IN.name)

# get matrix for spatialIDE
counts = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs_names)
coord = pd.DataFrame(adata.obsm['spatial'], columns=['x_coord', 'y_coord'], index=adata.obs_names)

# run spatialDE
results = SpatialDE.run(coord, counts)

# add results to .var
results.set_index("g", inplace=True)
results.index.name = None
results = results.add_prefix("spatialDE_")
adata.var = pd.concat([adata.var, results.loc[adata.var.index.values, :]], axis=1)

# write output
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
