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
    "output",
    type=argparse.FileType('w'),
    help='Output spatial h5ad file.'
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

# tangram specific modification
# add coordinates as obs.x and obs.y
if 'X_spatial' in adata_spatial.obsm:
    adata_spatial.obs['x'] = np.asarray(adata_spatial.obsm['X_spatial'][:,0])
    adata_spatial.obs['y'] = np.asarray(adata_spatial.obsm['X_spatial'][:,1])
elif 'spatial' in adata_spatial.obsm:
    adata_spatial.obs['x'] = np.asarray(adata_spatial.obsm['spatial'][:,0])
    adata_spatial.obs['y'] = np.asarray(adata_spatial.obsm['spatial'][:,1])
    adata_spatial.obsm['X_spatial'] = np.float32(adata_spatial.obsm['spatial'])[:,:2]
elif 'x' in adata_spatial.obs and 'y' in adata_spatial.obs:
    coord = pd.DataFrame({'x_coord': adata_spatial.obs['x'], 'y_coord': adata_spatial.obs['y']}, index=adata_spatial.obs_names)
    adata_spatial.obsm['X_spatial'] = np.float32(coord)
else:
    raise Exception("VSN ERROR: missing spatial coordinates in obsm or obs.")

# write output
adata_spatial.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
