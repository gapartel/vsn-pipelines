#!/usr/bin/env python3

import argparse
import numpy as np
import os
import scanpy as sc

parser = argparse.ArgumentParser(description='')

parser.add_argument(
    "input",
    nargs='+',
    type=argparse.FileType('r'),
    help='Input h5ad files.'
)

parser.add_argument(
    "-f", "--file-format",
    action="store",
    dest="format",
    default="h5ad",
    help="Concatenate the data. Choose one of : h5ad"
)

parser.add_argument(
    "-j", "--join",
    type=str,
    action="store",
    dest="join",
    default="inner",
    help="How to concatenate the multiple datasets. Choose one of : inner (intersect), outer (union)."
)

parser.add_argument(
    "-o", "--output",
    action="store",
    dest="output",
    default=None,
    help="Output file name."
)

args = parser.parse_args()

# Define the arguments properly
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output)[0]

# I/O
files = []
cell_ids = []

if args.format == 'h5ad':
    for FILE_PATH_IN in args.input:
        try:
            FILE_PATH_IN = FILE_PATH_IN.name
            adata = sc.read_h5ad(filename=FILE_PATH_IN)
            cell_ids.extend(adata.obs.index.values)
            files.append(adata)
        except IOError:
            raise Exception("VSN ERROR: Wrong input format. Expects .h5ad files, got .{}".format(FILE_PATH_IN))

index_unique = None

if len(cell_ids) != len(np.unique(cell_ids)):
    print("Non-unique cell index detected!")
    print("Make the index unique by joining the existing index names with the batch category, using index_unique='-'")
    index_unique = '-'
#
# Concatenate the data
#

if args.format == 'h5ad':
    # Concatenate multiple h5ad files
    # Source: https://anndata.readthedocs.io/en/latest/anndata.AnnData.concatenate.html#anndata.AnnData.concatenate
    adata = files[0].concatenate(
        files[1:],
        join=args.join,
        index_unique=index_unique
    )
    
    # If spatial dataset arrange 'X_spatial' embedding of each sample in a grid for SCope visualization
    if 'X_spatial' in adata.obsm.keys():
        samples = adata.obs.sample_id.unique()
        n_samples = len(samples)
        n_x = np.ceil(np.sqrt(n_samples)) if np.sqrt(n_samples) >= np.floor(np.sqrt(n_samples)) else np.floor(np.sqrt(n_samples))
        n_y = np.ceil(np.sqrt(n_samples)) if np.sqrt(n_samples) >= np.floor(np.sqrt(n_samples))+0.5 else np.floor(np.sqrt(n_samples))
       
        x = 0 
        x_max = y_max = 0
        X = []
        Y = []
        for s in samples:
            # Reset origin in (0,0)
            x_s = adata[adata.obs.sample_id==s].obsm['X_spatial'][:,0]-adata[adata.obs.sample_id==s].obsm['X_spatial'][:,0].min()
            y_s = adata[adata.obs.sample_id==s].obsm['X_spatial'][:,1]-adata[adata.obs.sample_id==s].obsm['X_spatial'][:,1].min()

            X.append(x_s + x_max)
            Y.append(y_s + y_max)

            x_max = x_max + x_s.max() + 1
            x = x + 1
            if x>n_x-1:
                x = 0
                x_max = 0
                y_max = y_max - y_s.max() - 1
        adata.obsm['X_spatial'] = np.array([np.concatenate(X), np.concatenate(Y)]).T

    # center data for SCope
    avg_coords = adata.obsm['X_spatial'].sum(axis=0)/adata.obsm['X_spatial'].shape[0]
    adata.obsm['X_spatial'] = adata.obsm['X_spatial'] - avg_coords
    # scale data for SCope such that X axis is between [-10,10]
    scfactor = 20/(np.max(adata.obsm['X_spatial'][:,0])-np.min(adata.obsm['X_spatial'][:,0]))
    adata.obsm['X_spatial'] = adata.obsm['X_spatial'] * scfactor

    # Not casting to float 64 bits can lead to not exact reproducible results. See:
    # - https://github.com/theislab/scanpy/issues/1612
    # - https://github.com/vib-singlecell-nf/vsn-pipelines/issues/295
    adata.X = adata.X.astype(np.float64)
    adata.var.index = adata.var.index.astype(str)
    adata = adata[:, np.sort(adata.var.index)]
    print(f"Total number of cells: {adata.obs.shape[0]}, genes: {adata.var.shape[0]}.")

   
else:
    raise Exception("VSN ERROR: Concatenation of .{} files is not implemented.".format(args.format))

# I/O
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))
