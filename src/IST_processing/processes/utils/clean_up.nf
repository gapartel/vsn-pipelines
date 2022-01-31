nextflow.enable.dsl=2

import java.nio.file.Paths

moduleName="utils"
binDir = Paths.get(workflow.projectDir.toString(), "src/IST_processing/bin/$moduleName/")


process clean_work_dir {

    script:
    """
    bash $binDir/clean_work.sh $params.global.outdir $workDir
    """
}
