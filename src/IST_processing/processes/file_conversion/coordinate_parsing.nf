nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName="file_conversion"

binDir = Paths.get(workflow.projectDir.toString(), "src/IST_processing/bin/$moduleName/")

process transform_tile_coordinate_system {
    publishDir "$params.global.outdir/final", mode: 'copy'

    input:
    path csv
    
    val grid_size_x
    val grid_size_y
    val tile_size_x
    val tile_size_y
    output:
    path "${csv.baseName}_transformed.csv"

    script:
    """
    python $binDir/transformTileCoordinateSystem.py $csv $grid_size_x $grid_size_y $tile_size_x $tile_size_y
    """
}

process transform_tile_coords_and_prune {
    publishDir "$params.global.outdir/outs", mode: 'copy'

    input:
    path csv
    
    val grid_size_x
    val grid_size_y
    val tile_size_x
    val tile_size_y
    output:
    path "coords.csv"

    script:
    """
    python $binDir/transformTileCoordinateSystem.py $csv $grid_size_x $grid_size_y $tile_size_x $tile_size_y prune
    """
}
