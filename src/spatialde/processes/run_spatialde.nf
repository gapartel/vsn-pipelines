nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/spatialde/bin/" : ""

process SC__SPATIALDE__VARIABLE_GENES {
    
    container params.sc.spatialde.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__cpu'

    input:
        tuple val(sampleId), path(f)

    output:
        tuple val(sampleId), path("${sampleId}.SC__SPATIALDE__VARIABLE_GENES.${processParams.off}")

    script:
	def sampleParams = params.parseConfig(sampleId, params.global, params.sc.spatialde)
      	    processParams = sampleParams.local
      	"""
     	${binDir}/run_spatialde.py \
		$f \
      		"${sampleId}.SC__SPATIALDE__VARIABLE_GENES.${processParams.off}"
      	"""
}
