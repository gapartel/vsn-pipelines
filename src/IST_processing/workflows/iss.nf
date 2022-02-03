nextflow.enable.dsl=2

///////////////////////
// Include processes:
///////////////////////

include {
    iss_round_adder;
} from './utils/image_name_parser.nf'
include {
    CLIP_AND_RESCALE_TILES; CLIP_AND_RESCALE_GLOBAL
} from "./normalization/normalization_workflow.nf"
include {
    standard_iss_tiling as tiling;
} from "./tiling/tiling_workflows.nf"
include {
    create_reference_image
} from "./utils/projections.nf"
include {
    white_tophat_filter
} from "./filtering/filter_workflow.nf"

include {
    local_registration_of_tiles as registering 
} from "./registration/local_registration.nf"

include {
    spot_detection_iss;
} from "./spot_detection/spot_detection.nf"

include {
    decode_sequential_max_intensity as decoding
} from "../processes/decoding/decoding.nf"

include {
    iss_decoding_statistics
} from "./analytics/decoded_statistics.nf"
include {
    assignment_statistics_workflow
} from "./analytics/assigned_statistics.nf"

include {
    plot_decoded_genes 
} from "./plotting/decoded_genes_workflow.nf" 
include {
   threshold_watershed_segmentation as segmentation //stardist_segmentation_workflow as segmentation //
} from "./segmentation/segmentation_workflow.nf"

include {
    transform_tile_coordinate_system
} from "../processes/file_conversion/coordinate_parsing.nf"

include {
    clean_work_dir
} from "../processes/utils/clean_up.nf"



workflow iss {
    take:
        rounds

    main:
        // Pipeline most definitely needs params.reference to exists, so create one out of the first round if it does not exist yet
        if (!params.data.containsKey("reference")){
            log.info "No Reference image found, one will be created by taking the maximum intensity projection of the first round."
            //If you ever want to remove the round tuple value from this:  rounds.groupTuple(by:0).map {round_nr, files -> files}.first()
            rounds = iss_round_adder()
            params.data.reference = create_reference_image(rounds.groupTuple(by:0).first()) //Create reference image by taking maxIP on the first round
        }
        // Create the channel of round images, reference is explicitely defined in the config file as params.reference
       rounds = Channel.fromPath("${params.data.dataDir}/${params.data.round_prefix}*/${params.data.round_prefix}*_${params.data.channel_prefix}*.${params.data.extension}")
       

       // Tiling workflow: including calculating the highest resolution, padded all images to a resolution that would tile all images in equal sizes, registering globally
       tiling("$params.data.dataDir/$params.data.round_prefix*/${params.data.round_prefix}*_${params.data.channel_prefix}*.$params.data.extension", rounds, params.data.reference, params.data.DAPI)
       // Output of this is used often, so we rename the global variables for readability:
       grid_size_x = tiling.out.grid_size_x
       grid_size_y = tiling.out.grid_size_y

       
       //perform white tophat filtering on both reference and round images
       white_tophat_filter(tiling.out.reference, tiling.out.rounds, grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)

       // Register tiles locally:
       registering(white_tophat_filter.out.filtered_ref, white_tophat_filter.out.filtered_round, grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)

        // Detect spots
       spot_detection_iss(white_tophat_filter.out.filtered_ref, registering.out, grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)
       
       // Decode spots
       decoding(spot_detection_iss.out)

       // Pool decoded genes into one file for downstream analysis
       decoding.out.collectFile(name: "$params.global.outdir/decoded/concat_decoded_genes.csv", sort:true, keepHeader:true).set {decoded_genes}

       transform_tile_coordinate_system(decoded_genes, grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y) // Add original X and Y coordinates for later downstream analysis

       // Plot decoded genes
       if (params.utils.plot==true){
           plot_decoded_genes(tiling.out.reference, decoding.out, decoded_genes, tiling.out.padded_whole_reference,  grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)
       }
       
       // Segment cells on dapi
       segmentation(tiling.out.dapi, decoding.out, tiling.out.reference,  grid_size_x, grid_size_y, params.data.target_tile_x, params.data.target_tile_y)

       // Get analytics from decoding
       iss_decoding_statistics(decoded_genes, decoding.out)

       assignment_statistics_workflow(segmentation.out.concat_assigned_genes)

       // Copy every symlink and clean up work dir 
       /* clean_work_dir() */
    
}
