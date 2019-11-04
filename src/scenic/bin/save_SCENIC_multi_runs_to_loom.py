#!/usr/bin/env python3

import argparse
import re
import sys
import pandas as pd
import pickle
from pyscenic import transform
from pyscenic.export import export2loom
from pyscenic.transform import COLUMN_NAME_NES
from pyscenic.utils import COLUMN_NAME_MOTIF_SIMILARITY_QVALUE, COLUMN_NAME_ORTHOLOGOUS_IDENTITY, \
    COLUMN_NAME_ANNOTATION
import time
import utils
import export_to_loom

################################################################################
# TODO:
# This implementation should be optimized:
# It's taking several hours to run (~5h for ~9k genes and ~13k cells)
################################################################################

parser_grn = argparse.ArgumentParser(description='Run AUCell on gene signatures saved as TSV in folder.')

parser_grn.add_argument(
    'expression_mtx_fname',
    type=argparse.FileType('r'),
    help='The name of the file that contains the expression matrix for the single cell experiment.'
         ' Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).'
)
parser_grn.add_argument(
    'motif_enrichment_table_fname',
    type=argparse.FileType('r'),
    help='The name of the file that contains the motif enrichments.'
)
parser_grn.add_argument(
    'signatures_fname',
    help='The name of the folder containing the signatures as TSV files.'
)
parser_grn.add_argument(
    'auc_mtx_fname',
    type=argparse.FileType('r'),
    help='The name of the file that contains the AUCell matrix.'
)
parser_grn.add_argument(
    '-o', '--output',
    type=argparse.FileType('w'),
    default=sys.stdout,
    help='Output file/stream, i.e. a table of TF-target genes (CSV).'
)
parser_grn.add_argument(
    '--min-genes-regulon',
    type=int,
    default=5,
    dest="min_genes_regulon",
    help='The threshold used for filtering the regulons based on the number of targets (default: {}).'.format(5)
)
parser_grn.add_argument(
    '--min-regulon-gene-occurrence',
    type=int,
    default=5,
    dest="min_regulon_gene_occurrence",
    help='The threshold used for filtering the genes bases on their occurrence (default: {}).'.format(5)
)
parser_grn.add_argument(
    '--cell-id-attribute',
    type=str,
    default='CellID',
    dest="cell_id_attribute",
    help='The name of the column attribute that specifies the identifiers of the cells in the loom file.'
)
parser_grn.add_argument(
    '--gene-attribute',
    type=str,
    default='Gene',
    dest="gene_attribute",
    help='The name of the row attribute that specifies the gene symbols in the loom file.'
)
parser_grn.add_argument(
    '--title',
    type=str,
    dest="title",
    help='The title for this loom file. If None than the basename of the filename is used as the title.'
)
parser_grn.add_argument(
    '--nomenclature',
    type=str,
    dest="nomenclature",
    help='The name of the genome.'
)
parser_grn.add_argument(
    '--scope-tree-level-1',
    type=str,
    dest="scope_tree_level_1",
    help='The name of the first level of the SCope tree.'
)
parser_grn.add_argument(
    '--scope-tree-level-2',
    type=str,
    dest="scope_tree_level_2",
    help='The name of the second level of the SCope tree.'
)
parser_grn.add_argument(
    '--scope-tree-level-3',
    type=str,
    dest="scope_tree_level_3",
    help='The name of the third level of the SCope tree.'
)

args = parser_grn.parse_args()

print(f"Extracting the matrix form the loom...", flush=True)
start = time.time()
ex_matrix_df = utils.get_matrix(
    loom_file_path=args.expression_mtx_fname.name,
    gene_attribute=args.gene_attribute,
    cell_id_attribute=args.cell_id_attribute
)
print(f"... took {time.time() - start} seconds", flush=True)

# Transform motif enrichment table (generated from the cisTarget step) to regulons
print(f"Reading aggregated motif enrichment table...", flush=True)
start = time.time()
f = args.motif_enrichment_table_fname.name
if f.endswith('.pickle'):
    with open(f, 'rb') as handle:
        motif_enrichment_table = pickle.load(handle)
elif f.endswith('.csv'):
    motif_enrichment_table = utils.read_feature_enrichment_table(fname=args.motif_enrichment_table_fname.name, sep=",")
else:
    raise Exception("The aggregated feature enrichment table is in the wrong format. Expecting .pickle or .csv formats.")
print(f"... took {time.time() - start} seconds to run.", flush=True)

print(f"Making the regulons...", flush=True)
start = time.time()
regulons = transform.df2regulons(
    df=motif_enrichment_table,
    save_columns=[
        COLUMN_NAME_NES,
        COLUMN_NAME_ORTHOLOGOUS_IDENTITY,
        COLUMN_NAME_MOTIF_SIMILARITY_QVALUE,
        COLUMN_NAME_ANNOTATION
    ]
)
print(f"{len(regulons)} regulons from df2regulons.")

# Read the signatures saved in out/multi_runs_regulons_[mtf|trk]
# Keep all regulons and targets (so that all can be visualized in SCope)
signatures = utils.read_signatures_from_tsv_dir(
    dpath=args.signatures_fname,
    noweights=False,
    weight_threshold=0,
    min_genes=0
)
print(f"{len(signatures)} all regulons from out/multi_runs_regulons_[mtf|trk].")

# Filter regulons (regulons from motifs enrichment table) by the filtered signatures
regulons = list(filter(lambda x: x.name in list(map(lambda x: x.name, signatures)), regulons))
# Add gene2occurrence from filtered signatures to regulons
regulons = list(
    map(
        lambda x:
        x.copy(gene2occurrence=list(filter(lambda y: y.name == x.name, signatures))[0].gene2weight), regulons
    )
)
print(f"{len(regulons)} final regulons.")
print(f"... took {time.time() - start} seconds to run.", flush=True)

print(f"Reading AUCell matrix...", flush=True)
start = time.time()
# Read the regulons AUCell matrix
auc_mtx = pd.read_csv(args.auc_mtx_fname.name, sep='\t', header=0, index_col=0)
auc_mtx.columns.name = "Regulon"
print(f"... took {time.time() - start} seconds to run.", flush=True)

# Create loom
print(f"Exporting to loom...", flush=True)
start = time.time()
# Create the basic loom
scope_loom = export_to_loom.SCopeLoom(
    ex_mtx=ex_matrix_df,
    regulons=regulons,
    out_fname=args.output.name,
    title=args.title,
    nomenclature=args.nomenclature,
    auc_mtx=auc_mtx,
    tree_structure=[args.scope_tree_level_1, args.scope_tree_level_2, args.scope_tree_level_3],
    compress=True,
    save_additional_regulon_meta_data=True
)
# Add additional stuff specific to multi-runs SCENIC
scope_loom.add_row_attr_regulon_gene_weights()
scope_loom.add_row_attr_regulon_gene_occurrences()
scope_loom.add_meta_data(_dict={
    "regulonSettings": {
        "min_genes_regulon": args.min_genes_regulon,
        "min_regulon_gene_occurrence": args.min_regulon_gene_occurrence
    }
})
scope_loom.export()

print(f"... took {time.time() - start} seconds to run.", flush=True)
print(f"Done.", flush=True)
