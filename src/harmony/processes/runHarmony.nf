nextflow.enable.dsl=2

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/harmony/bin/" : ""

process SC__HARMONY__HARMONY_MATRIX {
    
    container params.getToolParams("harmony").container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink'
    label 'compute_resources__default'

    input:
        tuple \
            val(sampleId), \
            path(f)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__HARMONY__HARMONY_MATRIX.tsv")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.getToolParams("harmony"))
		processParams = sampleParams.local
        varsUseAsArguments = processParams.varsUse.collect({ '--vars-use' + ' ' + it }).join(' ')
        """
        ${binDir}run_harmony.R \
            ${f} \
            --seed ${params.global.seed} \
            ${varsUseAsArguments} \
            --output-prefix "${sampleId}.SC__HARMONY__HARMONY_MATRIX"
        """

}
