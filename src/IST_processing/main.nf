nextflow.enable.dsl=2

include {
    add_parent_dir_to_file_name 
} from "./processes/utils/file_name_parsing.nf"

include{
    SPLIT_CZI_ROUNDS_INTO_CHANNEL_TIFS
} from "./workflows/file_conversion/czi_conversion.nf"

include {
    intensity_diagnosing
} from "./workflows/quality_control/intensity_workflows.nf"

include {
    iss as iss_pipeline
} from "./workflows/iss.nf"

include {
    merfish as merfish_pipeline
} from "./workflows/merfish.nf"


workflow rename_files{

    add_parent_dir_to_file_name()
}
workflow convert_czi {

    SPLIT_CZI_ROUNDS_INTO_CHANNEL_TIFS("$params.data.dataDir/*.czi")
}
workflow quality_control{

    intensity_diagnosing("$params.data.dataDir/$params.round_prefix*/${params.round_prefix}*_${params.channel_prefix}*.$params.data.extension")
}

workflow iss {
    take:
        data
    main:
        iss_pipeline(data)
    emit:
        iss_pipeline.out
        
}
workflow merfish {
    take:
        data
    main:
        merfish_pipeline(data)
}

