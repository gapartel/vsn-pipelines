nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName="normalization"
binDir = Paths.get(workflow.projectDir.toString(), "src/$moduleName/bin/")

process clip_and_rescale {
    publishDir "$params.outDir/normalized", mode: 'copy'

    input: 
    path image

    output:
    path "${image.baseName}_normalized.tif"

    script:
    """
    python $binDir/clipAndRescale.py $image $params.clip_percentile
    """
}
process equalize_histogram {
    publishDir "$params.outDir/normalized", mode: 'symlink'

    input: 
    path image

    output:
    path "${image.baseName}_equalized.tif"

    script:
    """
    python $binDir/equalizeHistogram.py $image
    """
}

process match_histogram {
    publishDir "$params.outDir/normalized", mode: 'symlink'

    input: 
    path reference
    path target

    output:
    path "${target.baseName}_matched.tif"

    script:
    """
    python $binDir/matchHistogram.py $reference $target
    """
}
