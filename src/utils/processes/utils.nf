nextflow.enable.dsl=2

import java.nio.file.Paths
import nextflow.config.ConfigParser
import static groovy.json.JsonOutput.*

binDir = !params.containsKey("test") ? "${workflow.projectDir}/src/utils/bin" : Paths.get(workflow.scriptFile.getParent().getParent().toString(), "utils/bin")

def boolean isCollectionOrArray(object) {    
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

def getToolParams(params, toolKey) {
    if(!toolKey.contains(".")) {
        return params[toolKey]
    }
    def entry = params
    toolKey.split('\\.').each { entry = entry?.get(it) }
    return entry
}

def isParamNull(param) {
    return param == null || param == "NULL"
}

def clean(params) {
   return params.findAll { !it.key.contains('-') }
}

def printMap(map) {
    prettyPrint(toJson(map))
}

def detectCellRangerVersionData = { cellRangerV2Data, cellRangerV3Data ->
    if(cellRangerV2Data.isDirectory() || cellRangerV3Data.isDirectory()) {
        if(cellRangerV2Data.exists()) {
            // Sanity checks
            if(new File(Paths.get(cellRangerV2Data.toString(), "genes.tsv.gz").toString()).exists())
                throw new Exception("VSN ERROR: Found genes.tsv.gz but expecting genes.tsv. The gene file should be uncompressed.")
            if(new File(Paths.get(cellRangerV2Data.toString(), "genes.tsv").toString()).exists())
                return [
                    version: 2,
                    path: cellRangerV2Data
                ]
            // Extract genome folder if a single one exists
            genomes = cellRangerV2Data.list()
            if(genomes.size() > 1 || genomes.size() == 0)
                throw new Exception("VSN ERROR: None or multiple genomes detected for the output generated by CellRanger v2. Selecting custom genome is currently not implemented.")
            genomeFilePath = Paths.get(cellRangerV2Data.toString(), genomes[0])
            // Sanity checks
            if(!new File(genomeFilePath.toString()).isDirectory())
                throw new Exception("VSN ERROR: Expecting a genome directory from the output generated by CellRanger v2.")
            if(new File(Paths.get(genomeFilePath.toString(), "genes.tsv.gz").toString()).exists())
                throw new Exception("VSN ERROR: Found compressed gene file (genes.tsv.gz) but expecting uncompressed gene file (genes.tsv). Use gunzip for instance to uncompress it.")
            if(new File(Paths.get(genomeFilePath.toString(), "barcodes.tsv.gz").toString()).exists())
                throw new Exception("VSN ERROR: Found compressed gene file (barcodes.tsv.gz) but expecting uncompressed gene file (barcodes.tsv). Use gunzip for instance to uncompress it.")
            if(new File(Paths.get(genomeFilePath.toString(), "matrix.mtx.gz").toString()).exists())
                throw new Exception("VSN ERROR: Found compressed gene file (matrix.mtx.gz) but expecting uncompressed gene file (matrix.mtx.gz). Use gunzip for instance to uncompress it.")
            if(!new File(Paths.get(genomeFilePath.toString(), "genes.tsv").toString()).exists())
                throw new Exception("VSN ERROR: Expecting a gene file genes.tsv file but none are found.")
            if(!new File(Paths.get(genomeFilePath.toString(), "barcodes.tsv").toString()).exists())
                throw new Exception("VSN ERROR: Expecting a barcode file barcodes.tsv file but none are found.")
            if(!new File(Paths.get(genomeFilePath.toString(), "matrix.mtx").toString()).exists())
                throw new Exception("VSN ERROR: Expecting a matrix file matrix.mtx file but none are found.")
            return [
                version: 2,
                path: file(genomeFilePath)
            ]
        } else if(cellRangerV3Data.exists()) {
            if(!new File(Paths.get(cellRangerV3Data.toString(), "features.tsv").toString()).exists() && !new File(Paths.get(cellRangerV3Data.toString(), "features.tsv.gz").toString()).exists())
                throw new Exception("VSN ERROR: Expecting either a features.tsv or features.tsv.gz file but none are found.")
            return [
                version: 3,
                path: cellRangerV3Data
            ]
        } else {
            throw new Exception("VSN ERROR: Cannot detect the version of the data format of CellRanger.")
        }
    } else {
        if(cellRangerV2Data.exists()) {
            return [
                version: 2,
                path: cellRangerV2Data
            ]
        } else if(cellRangerV3Data.exists()) {
            return [
                version: 3,
                path: cellRangerV3Data
            ]
        } else {
            throw new Exception("VSN ERROR: Cannot detect the version of the data format of CellRanger.")
        }
    }
}

def runPythonConverter = {
    processParams,
    sampleId,
    inputDataType,
    outputDataType,
    outputExtension,
    group,
    f ->
    return (
        """
        ${binDir}/sc_file_converter.py \
            --sample-id "${sampleId}" \
            ${!isParamNull(group) ? '--group-name ' + group.get(0) : ''} \
            ${!isParamNull(group) ? '--group-value ' + group.get(1) : ''} \
            ${processParams?.makeVarIndexUnique ? '--make-var-index-unique '+ processParams.makeVarIndexUnique : ''} \
            ${processParams?.tagCellWithSampleId ? '--tag-cell-with-sample-id '+ processParams.tagCellWithSampleId : ''} \
            ${processParams?.remove10xGEMWell ? '--remove-10x-gem-well '+ processParams.remove10xGEMWell : ''} \
            ${processParams?.useRaw ? '--use-raw '+ processParams.useRaw : ''} \
            --input-format $inputDataType \
            --output-format $outputDataType \
            ${f} \
            "${sampleId}.SC__FILE_CONVERTER.${outputExtension}"
        """
    )
}

def runRConverter = {
    processName,
    processParams,
    sampleId,
    inputDataType,
    outputDataType,
    outputExtension,
    group,
    f,
    sceMainLayer = null ->
    
    return (
        """
        ${binDir}/sc_file_converter.R \
            --sample-id "${sampleId}" \
            ${!isParamNull(group) ? '--group-name ' + group.get(0) : ''} \
            ${!isParamNull(group) ? '--group-value ' + group.get(1) : ''} \
            ${processParams?.tagCellWithSampleId ? '--tag-cell-with-sample-id '+ processParams.tagCellWithSampleId : ''} \
            ${processParams?.remove10xGEMWell ? '--remove-10x-gem-well '+ processParams.remove10xGEMWell : ''} \
            ${(processParams.containsKey('seuratAssay')) ? '--seurat-assay '+ processParams.seuratAssay : ''} \
            ${(processParams.containsKey('seuratMainLayer')) ? '--seurat-main-assay '+ processParams.seuratMainLayer : ''} \
            ${sceMainLayer != null ? '--sce-main-layer '+ sceMainLayer : ''} \
            --input-format $inputDataType \
            --output-format $outputDataType \
            ${f} \
            "${sampleId}.${processName}.${outputExtension}"
        """
    )
}

def getConverterContainer = { params, type ->
    switch(type) {
        case "cistopic":
            return params.tools.cistopic.container
        case "r":
            return "vibsinglecellnf/scconverter:0.0.1"
        break;
        case "python":
            return params.tools.scanpy.container
    }
}

def getConverter = { iff, off ->
    if(iff == "10x_atac_cellranger_mex_outs" && off == "cistopic_rds")
        return "cistopic"
    if((iff == "seurat_rds" && off == "h5ad")
        || (iff == "10x_cellranger_mex" && off == "sce_rds")
        || (iff == "sce_rds" && off == "h5ad"))
        return "r"
    return "python"
}

def getOutputExtension = { off ->
    switch(off) {
        case "cistopic_rds":
            return "cisTopic.Rds"
        case "sce_rds":
            return "SCE.Rds"
        case "h5ad":
            return "h5ad"
        default:
            throw new Exception("VSN ERROR: This output file format is not implemented.")
    }
    return "python"
}

process SC__FILE_CONVERTER {

    cache 'deep'
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    container "${getConverterContainer(params,converterToUse)}"
    label 'compute_resources__mem'

    input:
        tuple \
            val(sampleId), \
            path(f), \
            val(inputDataType), \
            val(outputDataType), \
            val(group)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__FILE_CONVERTER.${outputExtension}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.utils.file_converter)
        processParams = sampleParams.local

        switch(inputDataType) {
            case "10x_cellranger_mex":
                // Nothing to be done here
            break;
            case "10x_cellranger_mex_outs":
                // Reference: https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices-
                // Check if output was generated with CellRanger v2 or v3
                cellranger_outs_v2_mex = file("${f.toRealPath()}/${processParams.useFilteredMatrix ? "filtered" : "raw"}_gene_bc_matrices/")
                cellranger_outs_v3_mex = file("${f.toRealPath()}/${processParams.useFilteredMatrix ? "filtered" : "raw"}_feature_bc_matrix/")
                cellRangerData = detectCellRangerVersionData(cellranger_outs_v2_mex, cellranger_outs_v3_mex)
                f = cellRangerData.path
                inputDataType = "10x_cellranger_mex"
            break;
            case "10x_cellranger_h5":
                // Nothing to be done here
            break;
            case "10x_cellranger_h5_outs":
                // Check if output was generated with CellRanger v2 or v3
                cellranger_outs_v2_h5 = file("${f.toRealPath()}/${processParams.useFilteredMatrix ? "filtered" : "raw"}_gene_bc_matrices.h5")
                cellranger_outs_v3_h5 = file("${f.toRealPath()}/${processParams.useFilteredMatrix ? "filtered" : "raw"}_feature_bc_matrix.h5")
                cellRangerData = detectCellRangerVersionData(cellranger_outs_v2_h5, cellranger_outs_v3_h5)
                f = cellRangerData.path
                inputDataType = "10x_cellranger_h5"
            case "10x_atac_cellranger_mex_outs":
                // Nothing to be done here
            break;
            case "csv":
            case "tsv":
            case "h5ad":
            case "loom":
            case "seurat_rds":
                // Nothing to be done here
            break;
            
            default:
                throw new Exception("VSN ERROR: The given input format ${inputDataType} is not recognized.")
            break;
        }

        // Get the converter based on input file format and output file format
        converterToUse = getConverter(
            inputDataType,
            outputDataType
        )
        outputExtension = getOutputExtension(outputDataType)

        switch(converterToUse) {
            case "cistopic":
                """
                ${binDir}/create_cistopic_object.R \
                    --tenx_path ${f} \
                    --sampleId ${sampleId} \
                    --output ${sampleId}.SC__FILE_CONVERTER.${outputExtension}
                """
                break;
            case "r":
                runRConverter(
                    "SC__FILE_CONVERTER",
                    processParams,
                    sampleId,
                    inputDataType,
                    outputDataType,
                    outputExtension,
                    group,
                    f
                )
                break;
            case "python":
                runPythonConverter(
                    processParams,
                    sampleId,
                    inputDataType,
                    outputDataType,
                    outputExtension,
                    group,
                    f
                )
                break;
            default:
                throw new Exception("VSN ERROR: Unrecognized file converter.")
        }

}

process SC__FILE_CONVERTER_FROM_SCE {

    cache 'deep'
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    container "${getConverterContainer(params,converterToUse)}"
    label 'compute_resources__mem'

    input:
        tuple \
            val(sampleId), \
            path(f), \
            val(group)
        val(outputDataType)
        val(mainLayer)

    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.SC__FILE_CONVERTER_FROM_SCE.${outputDataType}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.utils.file_converter)
        processParams = sampleParams.local
        def _outputDataType = outputDataType
        converterToUse = getConverter(
            "sce_rds",
            _outputDataType
        )
        def outputExtension = getOutputExtension(_outputDataType)

        runRConverter(
            "SC__FILE_CONVERTER_FROM_SCE",
            processParams,
            sampleId,
            "sce_rds",
            _outputDataType,
            outputExtension,
            group,
            f,
            mainLayer
        )

}

process SC__FILE_CONCATENATOR {

    cache 'deep'
    container params.tools.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

    input:
        file("*")

    output:
        tuple val(params.global.project_name), path("${params.global.project_name}.SC__FILE_CONCATENATOR.${processParams.off}")

    script:
        processParams = params.utils.file_concatenator
        """
        ${binDir}/sc_file_concatenator.py \
            --file-format $processParams.off \
            ${(processParams.containsKey('join')) ? '--join ' + processParams.join : ''} \
            --output "${params.global.project_name}.SC__FILE_CONCATENATOR.${processParams.off}" *
        """

}

process SC__STAR_CONCATENATOR() {

    container params.tools.scanpy.container
    publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

    input:
        tuple \
            val(sampleId), \
            path(f)

    output:
        tuple \
            val(sampleId), \
            path("${params.global.project_name}.SC__STAR_CONCATENATOR.${processParams.off}")

    script:
        def sampleParams = params.parseConfig(sampleId, params.global, params.utils.star_concatenator)
        processParams = sampleParams.local
        id = params.global.project_name
        """
        ${binDir}/sc_star_concatenator.py \
            --stranded ${processParams.stranded} \
            --output "${params.global.project_name}.SC__STAR_CONCATENATOR.${processParams.off}" $f
        """

}

def getOutputFileName(params, tag, f, fileOutputSuffix, isParameterExplorationModeOn, stashedParams) {
    // - In parameter exploration mode,
    //   - returns 
    //     - <tag>.<fileOutputSuffix>.<stashedParams>.<f.extension> if fileOutputSuffix is null
    //     - <tag>.<stashedParams (joined)>.<f.extension> otherwise
    // - If not in parameter exporation mode AND
    //   - If pipelineOutputSuffix exists AND
    //   - If it is set to 'none'
    //     - returns:
    //       - <tag>.<f.extension>
    //       - <tag>.<pipelineOutputSuffix>.<f.extension>
    // - If fileOutputSuffix is null
    //   - returns:
    //     - <f.baseName>.<f.extension>
    // - In all other cases
    //   - returns:
    //     - <f.baseName>.<fileOutputSuffix>.<f.extension>
    if(isParameterExplorationModeOn && !isParamNull(stashedParams))
        return isParamNull(fileOutputSuffix) ? 
            "${tag}.${stashedParams.findAll { it != 'NULL' }.join('_')}.${f.extension}" :
            "${tag}.${fileOutputSuffix}.${stashedParams.findAll { it != 'NULL' }.join('_')}.${f.extension}"
    def utilsPublishParams = params.utils.publish
    if(utilsPublishParams?.pipelineOutputSuffix) {
        if(utilsPublishParams.pipelineOutputSuffix == 'none')
            return "${tag}.${f.extension}"
        if(utilsPublishParams.pipelineOutputSuffix.length() == 0)
            throw new Exception("VSN ERROR: The parameter 'params.utils.publish.outputFileSuffix' cannot be empty. If you don't want to add a suffix to the final output, please set this param to 'none'.")
        return utilsPublishParams.pipelineOutputSuffix
    }
    if(isParamNull(fileOutputSuffix))
        return "${f.baseName}.${f.extension}"
    return "${tag}.${fileOutputSuffix}.${f.extension}"
}

def getPublishDir = { outDir, toolName ->
    if(isParamNull(toolName))
        return "${outDir}/data"
    return "${outDir}/data/${toolName.toLowerCase()}"
}


process SC__PUBLISH {

    publishDir \
        "${getPublishDir(params.global.outdir,toolName)}", \
        mode: "${params.utils.publish?.mode ? params.utils.publish.mode: 'link'}", \
        saveAs: { filename -> "${outputFileName}" }

    label 'compute_resources__minimal'
    
    input:
        tuple \
            val(tag), \
            path(f), \
            val(stashedParams)
        val(fileOutputSuffix)
        val(toolName)
        val(isParameterExplorationModeOn)

    output:
        tuple \
            val(tag), \
            path(outputFileName), \
            val(stashedParams)

    script:
        outputFileName = getOutputFileName(
            params,
            tag,
            f,
            fileOutputSuffix,
            isParameterExplorationModeOn,
            stashedParams
        )
        /* avoid cases where the input and output files have identical names:
           Move the input file to a unique name, then create a link to
           the input file */
        """
        mv $f tmp
        if [ ! -f ${outputFileName} ]; then
            ln -L tmp "${outputFileName}"
        fi
        """
}


process COMPRESS_HDF5() {

	container "vibsinglecellnf/hdf5:1.10.5-r2"
	publishDir "${params.global.outdir}/data/intermediate", mode: 'symlink', overwrite: true
    label 'compute_resources__mem'

	input:
		tuple \
            val(id), \
            path(f), \
            val(stashedParams)
        val(fileOutputSuffix)
        val(toolName)

	output:
		tuple \
            val(id), \
            path("${outputFileName}"), \
            val(stashedParams)

	shell:
        def compressionLevel = params.utils?.publish &&
            params.utils.publish?.compressionLevel ? 
            params.utils.publish.compressionLevel : 
            6

        outputFileName = getOutputFileName(
            params,
            id,
            f,
            fileOutputSuffix,
            !isParamNull(stashedParams),
            stashedParams
        )
		"""
		GZIP_COMPRESSION_LEVEL=${compressionLevel}
        mv $f tmp
		h5repack \
		   -v \
		   -f GZIP=\${GZIP_COMPRESSION_LEVEL} \
		   tmp \
		   "${outputFileName}"
        rm tmp
		"""

}
