nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/tangram/bin/" : ""

process TANGRAM__MAP_CELLTYPES {
    
    container params.tools.tangram.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__gpu'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.TANGRAM__CELLTYPES.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.tools.tangram)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/run_tangram.py \
		-r  $processParams.reference_scrnaseq \
		${processParams?.device ? "-d " + processParams.device  : ""} \
		${processParams?.annotation ? "-a " + processParams.annotation  : ""} \
		${processParams?.gene_selection_method ? "-m " + processParams.gene_selection_method  : ""} \
		${processParams?.number_genes ? "-n " + processParams.number_genes  : ""} \
		${processParams?.qvalue ? "-q " + processParams.qvalue  : ""} \
		${processParams?.mapping_mode ? "--mode " + processParams.mapping_mode  : ""} \
		${processParams?.exp_scrnaseq ? "--exp-ref"  : ""} \
		${processParams?.exp_spatial ? "--exp-spatial"  : ""} \
		$f \
      		${sampleId}.TANGRAM__CELLTYPES.${processParams.off}
      	"""
}
