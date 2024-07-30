## Snakemake - bisantine-geo
##
## @alonzowolfram
##

# --- Necessary Python packages --- #
from datetime import datetime
import sys 

# --- Importing configuration files --- #
# https://stackoverflow.com/questions/67108673/accessing-the-path-of-the-configfile-within-snakefile
args = sys.argv
CONFIG_PATH = args[args.index("--configfiles") + 1]
configfile: CONFIG_PATH

# --- Setting variables --- #
# Output path
def generateOutputPath(previous_run_out_dir, output_path, project_name, run_name, now):
    # Set up the project_name and run_name strings.
    if project_name is None or project_name=="":
        project_name_string = ""
    else:
        project_name_string = project_name + "_"
    if run_name is None or run_name=="":
        run_name_string = ""
    else:
        run_name_string = run_name + "_"

    if previous_run_out_dir is None or previous_run_out_dir == "":
        if output_path is None or output_path == "":
            output_path = "out/" + project_name_string + run_name_string + str(now.year) + "-" + str(now.month) + "-" + str(now.day) + "-" + str(now.hour) + "-" + str(now.minute) + "-" + str(now.second) + "/"
        else:
            output_path = output_path
    else:
        output_path = previous_run_out_dir + "/"

    return output_path

now = datetime.now()
OUTPUT_PATH = generateOutputPath(config["data"]["previous_run_out_dir"], config["output"]["output_dir"], config["project"]["meta"]["project_name"], config["project"]["meta"]["run_name"], now)

# --- Rules --- # 
rule tcr_analysis: 
    input:
        script = "src/TCR_analysis.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds",
        previous_module = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_immune-deconvolution.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_TCR-analysis.rds",
        raw_plots = OUTPUT_PATH + "Rdata/TCR-analysis_plots-list.rds",
        anova_results = OUTPUT_PATH + "Rdata/TCR-analysis_ANOVA-res-list.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "tcr_analysis",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/TCR-analysis.out",
        err = OUTPUT_PATH + "logs/TCR-analysis.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule immune_deconvolution: 
    input:
        script = "src/immune_deconvolution.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds",
        previous_module = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_pathway-analysis.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_immune-deconvolution.rds",
        immune_deconv_results = OUTPUT_PATH + "Rdata/immune-deconvolution_results.rds",
        raw_plots = OUTPUT_PATH + "Rdata/immune-deconvolution_plots-list.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "immune_deconvolution",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/immune-deconvolution.out",
        err = OUTPUT_PATH + "logs/immune-deconvolution.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule pathway_analysis:
    input:
        script = "src/pathway_analysis.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_marker-identification.rds",
        DE_genes_table = OUTPUT_PATH + "tabular/LMM-differential-expression_results.csv"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_pathway-analysis.rds",
        pathways_table = OUTPUT_PATH + "tabular/pathway-analysis_results.csv"
    params:
        output_path = OUTPUT_PATH,
        current_module = "pathway_analysis",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/pathway-analysis.out",
        err = OUTPUT_PATH + "logs/pathway-analysis.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} {input.DE_genes_table} 1> {log.out} 2> {log.err}"

rule marker_identification:
    input:
        script = "src/marker_identification.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds",
        previous_module = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_differential-expression.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_marker-identification.rds",
        marker_table = OUTPUT_PATH + "tabular/LMM-marker_results.csv"
    params:
        output_path = OUTPUT_PATH,
        current_module = "marker_identification",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/marker-identification.out",
        err = OUTPUT_PATH + "logs/marker-identification.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule differential_expression_analysis:
    input:
        script = "src/differential_expression_analysis.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds",
        previous_module = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_unsupervised-analysis.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_differential-expression.rds",
        DE_genes_table = OUTPUT_PATH + "tabular/LMM-differential-expression_results.csv"
    params:
        output_path = OUTPUT_PATH,
        current_module = "differential_expression_analysis",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/differential-expression-analysis.out",
        err = OUTPUT_PATH + "logs/differential-expression-analysis.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule unsupervised_analysis:
    input:
        script = "src/unsupervised_analysis.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_unsupervised-analysis.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "unsupervised_analysis",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/unsupervised-analysis.out",
        err = OUTPUT_PATH + "logs/unsupervised-analysis.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule analysis_16s: 
    input:
        script = "src/16S_analysis.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_normalized.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "16S_analysis",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/16S-analysis.out",
        err = OUTPUT_PATH + "logs/16S-analysis.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule normalization:
    input:
        script = "src/normalization.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-probes.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_normalized.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "normalization",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/normalization.out",
        err = OUTPUT_PATH + "logs/normalization.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule qc_probes:
    input:
        script = "src/qc_probes.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-segments.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-probes.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "qc_probes",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/qc-probes.out",
        err = OUTPUT_PATH + "logs/qc-probes.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule qc_segments:
    input:
        script = "src/qc_segments.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-study-design.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-segments.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "qc_segments",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/qc-segments.out",
        err = OUTPUT_PATH + "logs/qc-segments.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule qc_study_design:
    input:
        script = "src/qc_study-design.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_raw.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-study-design.rds",
        Shiny_app = OUTPUT_PATH + "qc_segments_shiny_app.R"
    params:
        output_path = OUTPUT_PATH,
        current_module = "qc_study_design",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/qc-study-design.out",
        err = OUTPUT_PATH + "logs/qc-study-design.err" 
    shell:
        """
        Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}
        cp src/qc_segments_shiny_app.R {params.output_path}
        """

rule data_import_cleaning:
    input:
        script = "src/data_import_cleaning.R"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_raw.rds",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx"
    params:
        output_path = OUTPUT_PATH,
        current_module = "data_import_cleaning",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/data-import-cleaning.out",
        err = OUTPUT_PATH + "logs/data-import-cleaning.err" 
    shell:
        "Rscript {input.script} {params.config_path} {params.output_path} {params.current_module} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule export_env:
    output:
        env_file = OUTPUT_PATH + "config/conda.env"
    params:
        config_path = CONFIG_PATH,
        output_path = OUTPUT_PATH,
        config_file = CONFIG_PATH + "config.yaml"
    shell:
        "conda list --export > {params.output_path}config/conda.env; cp {params.config_path} {params.output_path}config/"