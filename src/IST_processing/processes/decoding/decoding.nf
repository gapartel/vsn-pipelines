nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName = "decoding"
binDir = Paths.get(workflow.projectDir.toString(), "src/IST_processing/bin/$moduleName/")


process decode_sequential_max_intensity {
    publishDir "$params.global.outdir/decoded", mode: 'symlink'

    input:
    path max_intensities

    output:
    path "decoded_tile*.csv"

    """
    python $binDir/decodeSequentialMaxIntensity.py ${max_intensities} ${params.data.codebook}
    """

}
process pixel_based_decoding {
    publishDir "$params.global.outdir/decoded", mode: 'symlink'

    input:
    val x_dim
    val y_dim
    tuple val(tile_nr), path(tile_images)

    output:
    path "decoded_${tile_nr}.csv"

    """
    python $binDir/decodePixelBased.py $x_dim $y_dim $tile_nr $params.data.codebook $params.data.bit_length $params.tools.IST_processing.distance_threshold $params.data.image_prefix $tile_images 
    """

}

process nn_pixel_based_decoding {
    publishDir "$params.global.outdir/decoded", mode: 'symlink'

    input:
    val x_dim
    val y_dim
    val min_area
    tuple val(tile_nr), path(tile_images)

    output:
    path "decoded_${tile_nr}.csv"

    """
    python $binDir/nnDecodePixelBased.py $x_dim $y_dim $min_area $tile_nr $params.data.codebook $params.data.bit_length $params.tools.IST_processing.distance_threshold $params.data.image_prefix $tile_images 
    """

}
