## Snakemake - bisantine-geo
##
## @alonzowolfram
##

# --- Importing configuration files --- #
configfile: "config.yaml"

# --- Necessary Python packages --- #
from datetime import datetime

# --- Setting variables --- #
def generateOutputPath(output_path, project_name, run_name, now):
    # Set up the project_name and run_name strings.
    if project_name is None or project_name=="":
        project_name_string = ""
    else:
        project_name_string = project_name + "_"
    if run_name is None or run_name=="":
        run_name_string = ""
    else:
        run_name_string = run_name + "_"

    if output_path is None or output_path == "":
        output_path = "out/" + project_name_string + run_name_string + str(now.year) + "-" + str(now.month) + "-" + str(now.day) + "-" + str(now.hour) + "-" + str(now.minute) + "-" + str(now.second) + "/"
    else:
        output_path = output_path

    return output_path

now = datetime.now()
OUTPUT_PATH = generateOutputPath(config["output"]["output_dir"], config["project"]["meta"]["project_name"], config["project"]["meta"]["run_name"], now)

# --- Rules --- # 
rule differential_expression_analysis:
    input:
        script = "src/differential_expression_analysis.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_unsupervised-analysis.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_differential-expression.rds",
        DE_genes_table = OUTPUT_PATH + "tabular/LMM-differential-expression_results.csv"
    params:
        output_path = OUTPUT_PATH,
        current_module = "differential_expression_analysis",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx"
    log:
        out = OUTPUT_PATH + "logs/differential_expression_analysis.out",
        err = OUTPUT_PATH + "logs/differential_expression_analysis.err" 
    shell:
        "Rscript {input.script} config.yaml {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule unsupervised_analysis:
    input:
        script = "src/unsupervised_analysis.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_normalized.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_unsupervised-analysis.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "unsupervised_analysis",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx"
    log:
        out = OUTPUT_PATH + "logs/unsupervised_analysis.out",
        err = OUTPUT_PATH + "logs/unsupervised_analysis.err" 
    shell:
        "Rscript {input.script} config.yaml {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule normalization:
    input:
        script = "src/normalization.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-probes.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_normalized.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "normalization",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx"
    log:
        out = OUTPUT_PATH + "logs/normalization.out",
        err = OUTPUT_PATH + "logs/normalization.err" 
    shell:
        "Rscript {input.script} config.yaml {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule qc_probes:
    input:
        script = "src/qc_probes.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-segments.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-probes.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "qc_probes",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx"
    log:
        out = OUTPUT_PATH + "logs/qc_probes.out",
        err = OUTPUT_PATH + "logs/qc_probes.err" 
    shell:
        "Rscript {input.script} config.yaml {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

rule qc_segments:
    input:
        script = "src/qc_segments.R",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_raw.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-segments.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "qc_segments",
        ppt_file = OUTPUT_PATH + "pubs/GeoMx-analysis_PowerPoint-report.pptx"
    log:
        out = OUTPUT_PATH + "logs/qc_segments.out",
        err = OUTPUT_PATH + "logs/qc_segments.err" 
    shell:
        "Rscript {input.script} config.yaml {params.output_path} {params.current_module} {input.R_file} {params.ppt_file} 1> {log.out} 2> {log.err}"

# rule qc_study_design:
#     input:
#         script = src/qc_study_design.R
#     output:
#         R_file = ,
#         ppt = 
#     shell:
#         "Rscript {input.script} config.yaml"

rule data_import_cleaning:
    input:
        script = "src/data_import_cleaning.R",
        env_file = OUTPUT_PATH + "config/conda.env",
        config_file = OUTPUT_PATH + "config/config.yaml"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_raw.rds"
    params:
        output_path = OUTPUT_PATH,
        current_module = "data_import_cleaning"
    log:
        out = OUTPUT_PATH + "logs/data_import_cleaning_stdout.out",
        err = OUTPUT_PATH + "logs/data_import_cleaning_stderr.err" 
    shell:
        "Rscript {input.script} config.yaml {params.output_path} {params.current_module} 1> {log.out} 2> {log.err}"

rule export_env:
    input: 
        config_file = "config.yaml"
    output:
        env_file = OUTPUT_PATH + "config/conda.env",
        config_file = OUTPUT_PATH + "config/config.yaml"
    params:
        output_path = OUTPUT_PATH
    shell:
        "conda list --export > {params.output_path}config/conda.env; cp config.yaml {params.output_path}config/"

# rule test:
#     input:
#         script = "src/setup.R"
#     params:
#         output_path = OUTPUT_PATH,
#         current_module = "test"
#     log:
#         out = OUTPUT_PATH + "logs/test_stdout.out",
#         err = OUTPUT_PATH + "logs/test_stderr.err" 
#     shell:
#         "Rscript {input.script} config.yaml {params.output_path} {params.current_module} 1> {log.out} 2> {log.err}"
