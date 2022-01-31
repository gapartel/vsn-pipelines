nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName= "analytics"
binDir = Paths.get(workflow.projectDir.toString(), "src/IST_processing/bin/$moduleName/")

/**
    Assignment analytics
**/
process get_assignment_stats {
    publishDir "$params.global.outdir/analytics/assigned/assets", mode: 'copy'

    input:
    path assigned_genes

    output:
    path "general_assignment_information.html"
    path "top10_assigned_cells.html"
    path "top10_assigned_genes.html"

    script:
    """
    python $binDir/extractStatsFromAssignedGenes.py $assigned_genes
    """
}

process create_assignment_html {
    publishDir "$params.global.outdir/analytics/assigned/", mode: 'copy'
    input:
    path template
    path general_assignment_information
    path top10_assigned_cells
    path top10_assigned_genes

    output: 
    path "assignment_report.html"

    script:
    """
    python $binDir/createAssignedHTMLreport.py $template $general_assignment_information $top10_assigned_cells $top10_assigned_genes
    """
}

/**
    Spot-based techniques analytics
**/
process get_spot_based_decoding_stats {
    publishDir "$params.global.outdir/analytics/decoded_stats/assets/", mode: 'copy'

    input:
    path decoded_genes

    output:
    path "general_stats.html"
    path "decoded_stat.html"
    path "recognized_barcodes_per_gene.html"
    path "unique_barcodes_called_counted.html"
    path "channels_called.html"
    path "barcodes_counted.png"
    path "tile_stats.html"
    path "recognized_genes_per_tile.png"
    /* env max_expressed_non_recognized_barcode, emit: most_prominent_unrecognized_barcode */

    script:

    """
    max_expressed_non_recognized_barcode=(`python $binDir/extractStatsFromDecodedBarcodes.py $decoded_genes $params.codebook $params.nr_rounds $params.nr_channels`)
    """
}

process plot_decoding_intensity_QC {
    publishDir "$params.global.outdir/analytics/decoded_stats/assets/", mode: 'copy'

    input:
    path decoded_genes

    output:
    path "decoding_intensity_QC.png"

    script:
    """
    python $binDir/plotDecodedIntensityQC.py $decoded_genes
    """
}

process create_spot_based_decoding_html {
    publishDir "$params.global.outdir/analytics/decoded_stats/", mode: 'copy'
    input:
    path template
    path general_stats
    path decoded_stat
    path recognized_barcodes_per_gene
    path unique_barcodes_called_counted
    path channels_called
    path barcodes_counted
    path tile_stats
    path recognized_genes_per_tile
    //decoding potential process
    path decoding_potential_plot
    // Decoding intensity qc
    path decoding_intensity_QC_plot

    output: 
    path "decoding_report.html"

    script:
    """
    python $binDir/createHTMLreport.py $template $general_stats $decoded_stat $recognized_barcodes_per_gene $unique_barcodes_called_counted $channels_called $barcodes_counted  $tile_stats $recognized_genes_per_tile $decoding_potential_plot $decoding_intensity_QC_plot
    """
}
/**
    pixel-based techniques analytics
**/
process get_pixel_based_decoding_stats {
    publishDir "$params.global.outdir/analytics/decoded_stats/assets/", mode: 'copy'

    input:
    path decoded_genes
    path codebook

    output:
    path "general_stats.html"
    path "top10_genes.html"
    path "bot10_genes.html"
    path "distributions.png"

    script:
    """
    python $binDir/extractStatsFromMERFISHDecodedBarcodes.py $decoded_genes $codebook
    """
}

process create_pixel_based_decoding_html {
    publishDir "$params.global.outdir/analytics/decoded_stats/", mode: 'copy'
    input:
    path template
    path general_stats
    path top10_genes
    path bot10_genes
    path distributions

    output: 
    path "decoding_report.html"

    script:
    """
    python $binDir/createDecodedMerfishHTMLreport.py $template $general_stats $top10_genes $bot10_genes $distributions
    """
}
