nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName="plotting"
binDir = Paths.get(workflow.projectDir.toString(), "src/$moduleName/bin/")

process plotDecodingPotential {
    publishDir "$params.outDir/analytics/decoded_stats/assets", mode: 'copy'
    input: 
    path decoded_genes

    output:
    path "decoding_potential_plot.png"

    script:
    """
    python $binDir/plotDecodingPotential.py $decoded_genes $params.codebook 
    """
}
process  plotTileDecodingPotential {
    publishDir "$params.outDir/analytics/decoded_stats/assets", mode: 'copy'
    input: 
    path decoded_genes

    output:
    path "${decoded_genes.baseName}_decoding_potential_plot.png"

    script:
    """
    python $binDir/plotTileDecodingPotential.py $decoded_genes $params.codebook 
    """

}
process plot_decoded_spots {
    publishDir "$params.outDir/plots", mode: 'copy'

    input:
    path decoded_genes
    path reference_image

    val grid_size_x
    val grid_size_y
    val tile_size_x
    val tile_size_y
    

    output:
    path "all_decoded_spots_plotted.png"
    path "decoded_genes_plotted.png"
    script:

    """
    python $binDir/plotDecodedGenes.py $reference_image $decoded_genes $grid_size_x,$grid_size_y $tile_size_x $tile_size_y
    """
}
process plot_detected_spots {
    publishDir "$params.outDir/plots", mode: 'copy'

    input:
    path detected_spots
    val grid_size_x
    val grid_size_y
    val tile_size_x
    val tile_size_y
    

    output:
    path "detected_spots_plotted.png"
    script:

    """
    python $binDir/plotDetectedSpots.py $detected_spots $grid_size_x,$grid_size_y $tile_size_x $tile_size_y
    """
}

process plot_detected_spots_on_tile {
    publishDir "$params.outDir/plots/tiles", mode: 'copy'

    input:
    tuple val(tile_nr), path(tile_image),path(detected_spots)
    
    
    output:
    path "${detected_spots.baseName}_plotted.png"
    script:
    """
    python $binDir/plotDetectedSpotsOnTile.py $tile_image $detected_spots 1
    """
}
process plot_decoded_genes_on_tile {
    publishDir "$params.outDir/plots/tiles", mode: 'copy'

    input:
    tuple val(tile_nr), path(tile_image), path(decoded_genes)
    
    output:
    path "${decoded_genes.baseName}_plotted.png"
    script:
    """
    python $binDir/plotDecodedGenesOnTile.py $tile_image $decoded_genes 1
    """
}

process plot_segmentation_labels_on_ref {
    publishDir "$params.outDir/plots/segmentation/", mode: 'copy'

    input:
    tuple val(tile_nr), path(labeled_image),path(original_image)
    output:

    path "${labeled_image.baseName}_overlay_REF.png"

    script:

    """
    python $binDir/plotLabeledImages.py $labeled_image $original_image REF
    """
}
process plot_segmentation_labels_on_dapi {
    publishDir "$params.outDir/plots/segmentation/", mode: 'copy'

    input:
    tuple val(tile_nr), path(labeled_image),path(original_image)
    output:
    path "${labeled_image.baseName}_overlay_DAPI.png"

    script:

    """
    python $binDir/plotLabeledImages.py $labeled_image $original_image DAPI 
    """
}

process plot_segmentation_labels {
    publishDir "$params.outDir/plots/segmentation/", mode: 'copy'

    input:
    path labeled_image
    output:
    path "${labeled_image.baseName}_plotted.png"

    script:

    """
    python $binDir/plotLabeledImages.py $labeled_image
    """
}

process plot_assigned_genes {
    publishDir "$params.outDir/plots/assigned_genes/", mode: 'copy'

    input:
    tuple val(tile_nr), path(assigned_genes), path(labeled_image)
    
    output:
    path "${assigned_genes.baseName}_plotted.png"

    script:
    """
    python $binDir/plotAssignedGenes.py $assigned_genes $labeled_image
    """
}

process plot_specific_barcode {
    publishDir "$params.outDir/plots/quality_control/", mode: 'copy'

    input:
    path image
    path decoded_genes
    val barcode
    
    
    output:
    path "${image.baseName}_{gene}_expression_plotted.png"

    script:
    """
    python $binDir/plotGeneExpression.py $image $decoded_genes $barcode
    """
}
