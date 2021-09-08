#!/usr/bin/env python3

import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Gene imputation from sn/scRNA-seq data.')

parser.add_argument(
    "spatial",
    type=argparse.FileType('r'),
    help='Input spatial h5ad file.'
)

parser.add_argument(
    "sc",
    type=argparse.FileType('r'),
    help='Input sn/scRNA-seq h5ad file.'
)

parser.add_argument(
    "knn_input",
    type=argparse.FileType('r'),
    help='knn-imputation result file.'
)

parser.add_argument(
    "genes",
    type=argparse.FileType('r'),
    help='csv file containin a list of genes to impute.'
)

parser.add_argument(
    "output",
    type=argparse.FileType('w'),
    help='Output spatial h5ad file with imputed genes.'
)


args = parser.parse_args()


# Define the arguments properly
FILE_PATH_SPATIAL = args.spatial
FILE_PATH_SC = args.sc
FILE_PATH_KNN_INPUT = args.knn_input
FILE_PATH_GENES = args.genes
FILE_PATH_OUT_BASENAME = os.path.splitext(args.output.name)[0]

# I/O
# Expects h5ad file
try:
    ad_spatial = sc.read_h5ad(filename=FILE_PATH_SPATIAL.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_SPATIAL)[0]))

try:
    ad_sc = sc.read_h5ad(filename=FILE_PATH_SC.name)
except IOError:
    raise Exception("Wrong input format. Expects .h5ad files, got .{}".format(os.path.splitext(FILE_PATH_SC)[0]))

try:
    knn_res = pickle.load( open( FILE_PATH_KNN_INPUT.name, "rb" ) )
except IOError:
    raise Exception("Wrong input format. Expects binary pickle files, got .{}".format(os.path.splitext(FILE_PATH_KNN_INPUT)[0]))

try:
    genes_df = pd.read_csv( FILE_PATH_GENES, names=['gene'], sep=',' )
except IOError:
    raise Exception("Wrong input format. Expects a single column .csv file with no header, got .{}".format(os.path.splitext(FILE_PATH_GENES)[0]))

#
# Run gene imputation
#
distances, indices = knn_res
genes_to_predict = genes_df.values
Imp_Genes = pd.DataFrame(np.zeros((ad_spatial.shape[0], genes_df.shape[0])),
                                 columns=genes_to_predict)

Spatial_data = ad_spatial.to_df()
RNA_data = ad_sc.to_df()

# Collect only those genes present in adata
genes_to_predict = genes_to_predict[np.intersect1d(RNA_data.columns,genes_to_predict)]
    
for j in range(0,Spatial_data.shape[0]):
    weights = 1-(distances[j,:][distances[j,:]<1])/(np.sum(distances[j,:][distances[j,:]<1]))
    weights = weights/(len(weights)-1)
    Imp_Genes.iloc[j,:] = np.dot(weights,RNA_data[genes_to_predict].iloc[indices[j,:][distances[j,:] < 1]])

# I/O
adata = ad_spatial.concatenate(
adata.write_h5ad("{}.h5ad".format(FILE_PATH_OUT_BASENAME))

