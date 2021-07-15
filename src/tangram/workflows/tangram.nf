nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:
// utils
include {
    PUBLISH as PUBLISH_H5AD_TANGRAM_CELLTYPES;
} from "../../utils/workflows/utils.nf" params(params)


include {
    TANGRAM__MAP_CELLTYPES;
} from '../processes/run_tangram.nf' params(params)

//////////////////////////////////////////////////////

workflow MAP_CELLTYPES {

    take:
        data

    main:
        TANGRAM__MAP_CELLTYPES( data )
        PUBLISH_H5AD_TANGRAM_CELLTYPES(
            TANGRAM__MAP_CELLTYPES.out.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "TANGRAM.celltypes_output",
            "h5ad",
            "tangram",
            false
        )

    emit:
        TANGRAM__MAP_CELLTYPES.out
}
