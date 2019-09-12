# SingleCellTxBenchmark

A repository of pipelines for single-cell data in Nextflow DSL2.

All the output generated by the pipelines will be located in the directory specified by the `--outdir` nextflow parameter.

# Dependencies

Make sure you have the following softwares installed,
- nextflow (version 19.08.1-edge.5132)

# Pipelines

## Multiple Datasets

Pipelines to aggregate multiple datasets together.

### bbknn 
Source: https://github.com/Teichlab/bbknn/blob/master/examples/pancreas.ipynb

**How to run on 10xGenomics datasets ?**

```{bash}
OUTPUT_DIRECTORY="out"
PROJECT_NAME="tiny"
```

Let's say the file structure of your data looks like this,

```
/home/data/
└── cellranger
    ├── Sample A
    │   └── outs
    │       ├── filtered_feature_bc_matrix
    │       └── ...
    └── Sample_B
        └── outs
            ├── filtered_feature_bc_matrix
            └── ...
```

Then the command to run the pipeline will be:

Using conda,
```
nextflow run \
   src/singlecelltxbenchmark/pipelines/bec__bbknn \
   -profile conda \
   --tenx_folder "/home/data/cellranger/**/filtered_feature_bc_matrix" \
   --sample_metadata /home/data/cellranger/metadata.tsv \
   --outdir ${OUTPUT_DIRECTORY} \
   --project_name ${PROJECT_NAME}
   --baseFilePath . \
   -with-report report.html \
   -with-trace
```

Using singularity,
```{bash}
nextflow run \
   src/singlecelltxbenchmark/pipelines/bec__bbknn \
      -profile singularity \
      --tenx_folder /home/data/cellranger/**/filtered_feature_bc_matrix \
      --sample_metadata /home/data/cellranger/metadata.tsv \
      --outdir ${OUTPUT_DIRECTORY} \
      --project_name ${PROJECT_NAME} \
      -with-report report.html \
      -with-trace
```