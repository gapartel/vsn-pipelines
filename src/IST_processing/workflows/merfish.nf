nextflow.enable.dsl=2

///////////////////////
// Include processes:
///////////////////////

include {
    rename_file
} from "../processes/utils/file_name_parsing.nf"
include {
    convert_to_uin16
} from "../processes/file_conversion/image_type_conversion.nf"

include {
    standard_merfish_tiling as tiling
} from "./tiling/tiling_workflows.nf"
include {
    merfish_global_registration
} from "./registration/merfish_registration.nf"
include {
    gaussian_low_pass_filter_workflow; gaussian_high_pass_filter_workflow; deconvolve_PSF_workflow; white_tophat_filter_merfish
} from "./filtering/filter_workflow.nf"
include {
    local_registration_of_tiles
} from "./registration/local_registration.nf"
include {
    merfish_threshold_watershed_segmentation as segmentation
} from "./segmentation/segmentation_workflow.nf"

include {
     nn_pixel_based_decoding as pixel_based_decoding  
} from "../processes/decoding/decoding.nf"

include {
    transform_tile_coordinate_system
} from "../processes/file_conversion/coordinate_parsing.nf"

include {
    plot_decoded_spots
} from "../processes/plotting/plotting.nf"

include{
    assignment_statistics_workflow
} from "./analytics/assigned_statistics.nf"
include {
        merfish_decoding_statistics 
} from "./analytics/decoded_statistics.nf"

workflow merfish {
    take:
        images

    main:
        // If n_tiles isn't in the config, that means the input images are whole-tissue images
        if (!params.data.containsKey("n_tiles")){
            glob_pattern ="${params.data.dataDir}/${params.data.image_prefix}*.${params.data.extension}" 
            images = convert_to_uin16(images)

            // Tile the images into equal sized tiles, and store the grid parameters
            // in variables for future use
            if (!params.containsKey("reference")){
                log.info "No Reference image found, one will be created by taking a random imaging round and renaming it."
                //If you even want to remove the round tuple value from this:  rounds.groupTuple(by:0).map {round_nr, files -> files}.first()
                params.data.reference = rename_file(images.first(), "REF") //Create reference image by taking maxIP on the first round
            }

            
            merfish_global_registration(params.data.reference, images)

            tiling(glob_pattern, merfish_global_registration.out, params.data.DAPI)
            grid_size_x = tiling.out.grid_size_x
            grid_size_y = tiling.out.grid_size_y
            tiled_rounds = tiling.out.rounds
            tiled_dapi = tiling.out.dapi
        }
        // Otherwise, the images are already tiled, so the input to image processing needs to be done differently
        else{
            glob_pattern ="${params.data.dataDir}/${params.data.tile_prefix}*${params.data.image_prefix}*.${params.data.extension}" 
            images = convert_to_uin16(images)
            tiled_rounds = images

            if (!params.data.containsKey("reference")){
                log.info "No Reference image found, one will be created by taking a random imaging round and renaming it."
                //If you even want to remove the round tuple value from this:  rounds.groupTuple(by:0).map {round_nr, files -> files}.first()
                params.data.reference = rename_file(images.first(), "REF") //Create reference image by taking maxIP on the first round
            }
            // We're assuming they are already registered and that they vollow the naming convention with tile*
            tiled_dapi =  Channel.fromPath(params.data.dapi_glob_pattern)
            grid_size_x = params.data.grid_size_x
            grid_size_y = params.data.grid_size_y
        }
        
        
        // Gaussian high pass filter 
        gaussian_high_pass_filter_workflow(tiled_rounds, grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)

        // Deconvolution with richardson_lucy
        deconvolve_PSF_workflow(gaussian_high_pass_filter_workflow.out, grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)

        gaussian_low_pass_filter_workflow(deconvolve_PSF_workflow.out, grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)
        // Map the images to their respective tiles, since for decoding they need to be in the correct order
        gaussian_low_pass_filter_workflow.out.map {file -> tuple((file.baseName=~ /tile\d+/)[0], file)} \
                                        | groupTuple()
                                        | set {grouped_images}
                                        
        pixel_based_decoding(params.data.target_tile_x, params.data.target_tile_y, params.tools.IST_processing.min_area, grouped_images)
        pixel_based_decoding.out.collectFile(name: "$params.global.outdir/decoded/concat_decoded_genes.csv", sort:true, keepHeader:true).set {decoded_genes}

        // Analysis processes only to do if the input image are whole-slide
        if (!params.data.containsKey("n_tiles")){
            transform_tile_coordinate_system(decoded_genes, grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)
            plot_decoded_spots(decoded_genes, tiling.out.padded_whole_reference, grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)
        }

        merfish_decoding_statistics(decoded_genes, params.data.codebook)

        segmentation(tiled_dapi, pixel_based_decoding.out, grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)
        // Calculate assignment stats
        assignment_statistics_workflow(segmentation.out.concat_assigned_genes)

        emit:
             segmentation.out.count_matrix


}
