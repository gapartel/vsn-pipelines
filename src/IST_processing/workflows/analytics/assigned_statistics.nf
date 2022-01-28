:nextflow.enable.dsl=2

include{
        get_assignment_stats ;  create_assignment_html 
} from "../../processes/analytics/analytics.nf"

workflow assignment_statistics_workflow{
        take:
            assigned_genes
        main:
            // General statistics
            get_assignment_stats(assigned_genes)
            
            create_assignment_html("$projectDir/src/IST-processing/assets/html_templates/assignment_report_template.html", get_assigned_stats.out)
}
