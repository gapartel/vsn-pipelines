nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName="file_conversion"

binDir = Paths.get(workflow.projectDir.toString(), "src/IST_processing/bin/$moduleName/")

process convert_to_uin16 {
    publishDir "$params.global.outdir/converted", mode: 'copy'

    input:
    path image
    
    output:
    path "${image.baseName}_16bit.tif"

    script:
    """
    python $binDir/convertToUint16.py $image
    """
}
