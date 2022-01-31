nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName="utils"
binDir = Paths.get(workflow.projectDir.toString(), "src/IST_processing/bin/$moduleName/")

process maxIP {
    publishDir "$params.global.outdir/maxIP/", mode: 'symlink'
    echo = true
    input:
    tuple val(round_nr), path(images)

    output:
    path "${round_nr}_maxIP.tif"


    script:
    """
    python $binDir/maxIP.py $round_nr $images
    """
}
