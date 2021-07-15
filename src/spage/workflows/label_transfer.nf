nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//  process imports:
// utils
include {
    PUBLISH as PUBLISH_H5AD_LABEL_TRANSFER;
} from "../../utils/workflows/utils.nf" params(params)

include {
    DATA_INTEGRATION;
} from '../processes/data_integration.nf' params(params)

include {
    KNN_IMPUTATION;
} from '../processes/knn_imputation.nf' params(params)

include {
    LABEL_TRANSFER;
} from '../processes/label_transfer.nf' params(params)

//////////////////////////////////////////////////////

workflow SPAGE__LABEL_TRANSFER {

    take:
        spatial_data
        sc_data
    main:
        data_integration_res = DATA_INTEGRATION( spatial_data, sc_data )
        knn_imputation_res = KNN_IMPUTATION( data_integration_res )
        label_transfer_res = LABEL_TRANSFER( data_integration_res.join(knn_imputation_res) )


        PUBLISH_H5AD_LABEL_TRANSFER(
            label_transfer_res.map {
                // if stashedParams not there, just put null 3rd arg
                it -> tuple(it[0], it[1], it.size() > 2 ? it[2]: null)
            },
            "SPAGE.label_transfer_output",
            "h5ad",
            "spage",
            false
        )
    emit:
        scanpyh5ad = label_transfer_res
}
