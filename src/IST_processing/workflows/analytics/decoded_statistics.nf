nextflow.enable.dsl=2


include{
        get_spot_based_decoding_stats  ; create_spot_based_decoding_html ; plot_decoding_intensity_QC
} from "../../processes/analytics/analytics.nf"

include {
       get_pixel_based_decoding_stats; create_pixel_based_decoding_html
} from "../../processes/analytics/analytics.nf"

include {
    plotDecodingPotential; plotTileDecodingPotential
} from "../../processes/plotting/plotting.nf"

workflow iss_decoding_statistics{
        take:
            decoded_genes
            decoded_genes_per_tile
        main:
            // General statistics
            get_spot_based_decoding_stats(decoded_genes)

            // Decoding potential throughout round progression
            plotDecodingPotential(decoded_genes)


            // Decoding intensity based on thresholds
            plot_decoding_intensity_QC(decoded_genes)
            
            create_spot_based_decoding_html("$projectDir/src/IST-processing/assets/html_templates/decoding_report_template.html",get_spot_based_decoding_stats.out, plotDecodingPotential.out, plot_decoding_intensity_QC.out)

}


workflow merfish_decoding_statistics{
        take:
            decoded_genes
            codebook

        main:
            // General statistics
            get_pixel_based_decoding_stats(decoded_genes, codebook)

            create_pixel_based_decoding_html("$projectDir/src/IST-processing/assets/html_templates/merfish_decoding_report_template.html",get_pixel_based_decoding_stats.out)

}
