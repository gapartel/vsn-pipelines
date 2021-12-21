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
    help='Input spatial h5ad file.'
)

parser.add_argument(
    "filtered",
    type=argparse.FileType('r'),
    help='Input raw filtered spatial h5ad file.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output spatial h5ad file.'
)


args = parser.parse_args()

# Define the arguments properly
FILE_PATH_IN = args.input
FILE_PATH_FILTERED = args.filtered
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

### main

# I/O
# Expects h5ad file
try:
    adata = sc.read_h5ad(filename=FILE_PATH_IN.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")

# Expects h5ad file
try:
    adata_raw = sc.read_h5ad(filename=FILE_PATH_FILTERED.name)
except IOError:
    raise Exception("VSN ERROR: Can only handle .h5ad files.")


# tangram specific modification
# add coordinates as obs.x and obs.y
adata_spatial = adata.copy()
adata_spatial.obs['x'] = np.asarray(adata.obsm['spatial'][:,0])
adata_spatial.obs['y'] = np.asarray(adata.obsm['spatial'][:,1])

# add raw spatial data
adata_spatial.X = adata_raw.X.copy()

# write output
adata_spatial.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
