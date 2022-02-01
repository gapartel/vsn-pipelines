nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName="filtering"
binDir = Paths.get(workflow.projectDir.toString(), "src/IST_processing/bin/$moduleName/")


process filter_ref {
    publishDir "$params.global.outdir/filtered_ref/", mode: 'symlink'

    input:
    path image 
    output:
    path "${image.baseName}_filtered.tif" 

    """
    python $binDir/whiteTophatFilter.py ${image} ${params.tools.IST_processing.filter_radius}
    """
}

process filter_round{
    // echo true
    publishDir "$params.global.outdir/filtered_round/", mode: 'symlink'
    
    input: 
    path image 

    output:
    path "${image.baseName}_filtered.tif"

    script:
    """
    python $binDir/whiteTophatFilter.py ${image} ${params.tools.IST_processing.filter_radius}
    """
}
process filter_gaussian_high_pass{
    publishDir "$params.global.outdir/filtered/high_passed/", mode: 'symlink'
    
    input: 
    path image 

    output:
    path "${image.baseName}_high_passed.tif"

    script:
    """
    python $binDir/gaussianHighPass.py ${image} ${params.tools.IST_processing.high_pass_sigma}
    """
}

process filter_gaussian_low_pass{
    publishDir "$params.global.outdir/filtered/low_passed", mode: 'symlink'
    
    input: 
    path image 

    output:
    path "${image.baseName}_low_passed.tif"

    script:

    """
    python $binDir/gaussianLowPass.py ${image} ${params.tools.IST_processing.low_pass_sigma}
    """
}

process deconvolve_PSF {
    publishDir "$params.global.outdir/filtered/deconvolved/", mode: 'symlink'
    
    input: 
    path image 

    output:
    path "${image.baseName}_deconvolved.tif"

    script:
    """
    python $binDir/restoreRichardsonLucy.py ${image} $params.tools.IST_processing.deconvolve_sigma $params.tools.IST_processing.iterations 
    """
}
