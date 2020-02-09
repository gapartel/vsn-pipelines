nextflow.preview.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/scanpy/bin/" : ""

import groovy.transform.TupleConstructor
import groovyx.gpars.dataflow.DataflowBroadcast

include '../../utils/processes/utils.nf'

@TupleConstructor()
class SC__SCANPY__CLUSTERING_PARAMS {

	Script env = null;
	LinkedHashMap configParams = null;
	// Parameters definiton
	String iff = null;
	String off = null;
	String report_ipynb = null;
	// Parameters benchmarkable
	String clusteringMethod = null; ArrayList<String> clusteringMethods = null;
    Float resolution = null; ArrayList<Float> resolutions = null;

	void setEnv(env) {
		this.env = env
	}

	void setConfigProcessParams(params) {
		this.configProcessParams = params
	}

	void displayMessage(tag) {
		Channel.from('').view {
			"""
------------------------------------------------------------------
\u001B[32m Benchmarking SC__SCANPY__CLUSTERING step... \u001B[0m
\u001B[32m Tag: ${tag} \u001B[0m
\u001B[32m Parameters tested: \u001B[0m
\u001B[32m - method: \u001B[0m \u001B[33m     ${clusteringMethods instanceof List} \u001B[0m
\u001B[32m   - values: \u001B[0m \u001B[33m   ${clusteringMethods} \u001B[0m
\u001B[32m - resolution: \u001B[0m \u001B[33m ${resolutions instanceof List} \u001B[0m
\u001B[32m   - values: \u001B[0m \u001B[33m   ${resolutions} \u001B[0m
------------------------------------------------------------------
            """
        }
	}

	String getMethodAsArgument(method) {
		return !this.env.isParamNull(method) ? '--method ' + method : ''
	}

	String getResolutionAsArgument(resolution) {
		return !this.env.isParamNull(resolution) ? '--resolution ' + resolution : ''
	}

	int numParamsBenchmarked() {
		def paramsBenchmarked = [ 
			clusteringMethods instanceof List,
			resolutions instanceof List
		]
		def sum = { result, i -> result + (i ? 1 : 0) }
		return paramsBenchmarked.inject(0, sum)
	}

	int numParams() {
		return 2 // Total number of parameters implemented for benchmarking
	}

	// Define a function to check if the current process is running in benchmark mode
	boolean isBenchmarkMode() {
		return (clusteringMethods instanceof List
				|| resolutions instanceof List
		)
	}

	DataflowBroadcast $(tag) {
		// Prepare argument stream
		def $method = Channel.from(clusteringMethods == null ? "NULL" : clusteringMethods)
		def $resolution = Channel.from(resolutions == null ? "NULL" : resolutions)
		displayMessage(tag)
		return $method.combine($resolution)
	}

}

def SC__SCANPY__CLUSTERING_PARAMS(params) {
	return (new SC__SCANPY__CLUSTERING_PARAMS(params))
}

/**
 * DEFAULT VERSION OF SCANPY CLUSTERING
 */
process SC__SCANPY__CLUSTERING {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true

  	input:
    	tuple val(sampleId), path(f)

  	output:
    	tuple val(sampleId), path("${sampleId}.SC__SCANPY__CLUSTERING.${processParams.off}")

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.clustering)
		processParams = sampleParams.local
		"""
		${binDir}cluster/sc_clustering.py \
			${(processParams.containsKey('clusteringMethod')) ? '--method ' + processParams.clusteringMethod : ''} \
			${(processParams.containsKey('resolution')) ? '--resolution ' + processParams.resolution : ''} \
			$f \
			"${sampleId}.SC__SCANPY__CLUSTERING.${processParams.off}"
		"""

}

/**
 * BENCHMARK VERSION OF SCANPY CLUSTERING
 */
process SC__SCANPY__BENCHMARK_CLUSTERING {

  	container params.sc.scanpy.container
  	clusterOptions "-l nodes=1:ppn=2 -l pmem=30gb -l walltime=1:00:00 -A ${params.global.qsubaccount}"
  	publishDir "${params.global.outdir}/data/intermediate/clustering/${isParamNull(method) ? "default": method.toLowerCase()}/${isParamNull(resolution) ? "default" : "res_" + resolution}", mode: 'symlink', overwrite: true

  	input:
    	tuple \
			val(sampleId), \
			path(f), \
			val(inertParams), \
			val(method), \
			val(resolution)

  	output:
    	tuple \
			val(sampleId), \
			path("${sampleId}.SC__SCANPY__BENCHMARK_CLUSTERING.${processParams.off}"), \
			val(method), \
			val(resolution)

  	script:
		def sampleParams = params.parseConfig(sampleId, params.global, params.sc.scanpy.clustering)
		processParams = sampleParams.local
		def _processParams = new SC__SCANPY__CLUSTERING_PARAMS()
		_processParams.setEnv(this)
		_processParams.setConfigParams(processParams)
		"""
		${binDir}cluster/sc_clustering.py \
			${_processParams.getMethodAsArgument(method)} \
			${_processParams.getResolutionAsArgument(resolution)} \
			$f \
			"${sampleId}.SC__SCANPY__BENCHMARK_CLUSTERING.${processParams.off}"
		"""

}
