nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:
// utils
include {
    PUBLISH as PUBLISH_H5AD_SPATIAL_VARIABLE_GENES;
} from "../../utils/workflows/utils.nf" params(params)


include {
    SC__SPATIALDE__VARIABLE_GENES;
} from '../processes/run_spatialde.nf' params(params)

//////////////////////////////////////////////////////

workflow GET_SPATIAL_VARIABLE_GENES {

    take:
        data

    main:
        SC__SPATIALDE__VARIABLE_GENES( data )
        PUBLISH_H5AD_SPATIAL_VARIABLE_GENES(
            SC__SPATIALDE__VARIABLE_GENES.out.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "SPATIALDE.variable_genes_output",
            "h5ad",
            "spatialDE",
            false
        )

    emit:
        SC__SPATIALDE__VARIABLE_GENES.out
}
