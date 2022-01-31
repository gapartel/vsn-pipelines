nextflow.enable.dsl=2

include {
    spot_detection_reference;spot_detection_round ; gather_intensities_in_rounds; get_max_intensities_over_channels
} from "../../processes/spot_detection/spot_detection.nf"

include {
    calculate_iss_precision_and_recall 
} from "../../workflows/quality_control/spot_detection_qc_workflows.nf"

include {
    transform_tile_coordinate_system ; transform_tile_coordinate_system as transform_tile_coordinate_system2
} from "../../processes/file_conversion/coordinate_parsing.nf"
include {
    plot_spots_whole_and_on_tiles
} from "../../workflows/plotting/detected_spots_workflow.nf" 


workflow spot_detection_iss {
    take:
    // Data
    reference //Preferably filtered by a filtering method: spot detection will be performed on this
    round_images //the detected spot's intensity will be measured on these images

    // Tile grid parameter
    grid_size_x
    grid_size_y
    tile_size_x
    tile_size_y

    main:
    //detect spots on reference image
    spot_detection_reference(reference)

    ////////////////////////////////////////////////////////
    //This is for spot detection quality control purposes //
    ////////////////////////////////////////////////////////
    if (params.tools.IST_processing.spot_detectionQC==true){
        transform_tile_coordinate_system(spot_detection_reference.out, grid_size_x, grid_size_y, tile_size_x, tile_size_y).set{transformed_ref_spots}
        transform_tile_coordinate_system.out.collectFile(name: "$params.global.outdir/blobs/transformed_concat_blobs.csv", sort:true, keepHeader:true).set {transformed_blobs}

        spot_detection_round(round_images)
        spot_detection_round.out.collectFile(name: "$params.global.outdir/hybs/concat_hybs.csv", sort:true, keepHeader:true).set {hybs}

        transform_tile_coordinate_system2(spot_detection_round.out, grid_size_x, grid_size_y, tile_size_x, tile_size_y) .set {transformed_round_spots}

        calculate_iss_precision_and_recall(transformed_blobs, transformed_round_spots) 
    }
    ///////////// end spot detection QC ///////////////////


    // Collect all spots in a seperate file
    spot_detection_reference.out.collectFile(name: "$params.global.outdir/blobs/concat_blobs.csv", sort:true, keepHeader:true).set {blobs}
    blobs_value_channel = blobs.first() //Call first to convert it into a value channel to allow for multiple iteration of a process with multiple inputs

    if (params.utils.plot==true){
        // Plot all detected spots
        plot_spots_whole_and_on_tiles(spot_detection_reference.out, blobs, reference, grid_size_x,grid_size_y, tile_size_x, tile_size_y)
    }
    
    //map round images into a tuple containing their round, channel and tile number, gather intensity code requires this information up front
    round_images.map {file -> tuple((file.baseName=~ /tile\d+/)[0],(file.baseName=~ /Round\d+/)[0],(file.baseName=~ /c\d+/)[0], file) }.set {round_images_mapped}

    // Gather intesnity on each round image of every spot
    gather_intensities_in_rounds(blobs_value_channel, round_images_mapped)
    // Collect all data into one file for downstream analysis
    gather_intensities_in_rounds.out.collectFile(name: "$params.global.outdir/intensities/concat_intensities.csv", sort:true, keepHeader:true).set {intensities}
    intensities_value_channel = intensities.first()

    get_max_intensities_over_channels(intensities_value_channel) 
    


    emit:
        get_max_intensities_over_channels.out.flatten()


}
