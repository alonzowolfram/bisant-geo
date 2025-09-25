## Snakemake - bisantine-geo
##
## @alonzowolfram
##

# --- Necessary Python packages --- #
from datetime import datetime
import sys 
import os
import filecmp
import shutil

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

# Nextflow vs Snakemake stuff
WORKFLOW_SYSTEM = "Snakemake"
PROJECT_DIRECTORY = ""

# --- Rules --- # 
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#onstart-onsuccess-and-onerror-handlers
onsuccess:
    workflow_system = WORKFLOW_SYSTEM,
    config_path = CONFIG_PATH,
    output_path = OUTPUT_PATH,
    project_directory = PROJECT_DIRECTORY,
    out = output_path[0] + "logs/make_report.out",
    err = output_path[0] + "logs/make_report.err",
    R_file = output_path[0] + "Rdata/latest_rule.Rds",
        
    # Ensure the output directory exists.
    os.makedirs(output_path[0] + "config/", exist_ok=True)

    # Get the current timestamp.
    timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

    # Export the current conda environment.
    os.system("conda list --export > " + output_path[0] + "config/conda_" + timestamp + ".env")

    # Extract the base name of the input config file (without extension.)
    config_name, config_ext = os.path.splitext(os.path.basename(config_path[0]))

    # Path to the new file.
    new_file = os.path.join(output_path[0] + "config/", f"{config_name}_{timestamp}{config_ext}")

    # Find the latest file in the `config` directory that matches this config base name
    existing_files = sorted(
        [f for f in os.listdir(output_path[0] + "config/") if f.startswith(config_name) and f.endswith(config_ext)]
    )
    latest_file = os.path.join(output_path[0] + "config/", existing_files[-1]) if existing_files else None

    # Compare the current YAML file with the latest file in the folder.
    if not latest_file or not filecmp.cmp(config_path[0], latest_file, shallow=False):
        # Copy the file if there are differences or no previous file exists.
        shutil.copy(config_path[0], new_file)
        print(f"Copied {config_path[0]} to {new_file}")
    else:
        print("No changes detected in configuration YAML file. Skipping copy.")

    # Create the report. 
    os.system("Rscript src/make_report.R " + config_path[0] + " " + workflow_system[0] + " " + "make_report" + " " + output_path[0] + " " + R_file[0] + " " + project_directory[0] + " 1> " + out[0] + " 2> " + err[0])

rule tcr_analysis: 
    input:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds",
        previous_module = OUTPUT_PATH + "Rdata/immune-deconvolution_results.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_TCR-analysis.rds"
        # ,raw_plots = OUTPUT_PATH + "Rdata/TCR-analysis_plots-list.rds",
        # anova_results = OUTPUT_PATH + "Rdata/TCR-analysis_ANOVA-res-list.rds"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/TCR_analysis.R",
        output_path = OUTPUT_PATH,
        current_module = "TCR_analysis",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/TCR-analysis.out",
        err = OUTPUT_PATH + "logs/TCR-analysis.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} 1> {log.out} 2> {log.err}"

rule immune_deconvolution: 
    input:
        previous_module = OUTPUT_PATH + "tabular/pathway-analysis_results.csv",
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds"
    output:
        immune_deconv_results = OUTPUT_PATH + "Rdata/immune-deconvolution_results.rds",
        raw_plots = OUTPUT_PATH + "Rdata/immune-deconvolution_plots-list.rds"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/immune_deconvolution.R",
        output_path = OUTPUT_PATH,
        current_module = "immune_deconvolution",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/immune-deconvolution.out",
        err = OUTPUT_PATH + "logs/immune-deconvolution.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} 1> {log.out} 2> {log.err}"

rule pathway_analysis:
    input:
        DE_genes_table = OUTPUT_PATH + "tabular/LMM-differential-expression_results.csv"
    output:
        pathways_table = OUTPUT_PATH + "tabular/pathway-analysis_results.csv"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/pathway_analysis.R",
        output_path = OUTPUT_PATH,
        current_module = "pathway_analysis",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/pathway-analysis.out",
        err = OUTPUT_PATH + "logs/pathway-analysis.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.DE_genes_table} 1> {log.out} 2> {log.err}"

rule differential_expression_analysis:
    input:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds",
        previous_module = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_unsupervised-analysis.rds"
    output:
        R_file_differential_expression_plot_list = OUTPUT_PATH + "Rdata/LMM-DEG_volcano-plots.rds",
        R_file_differential_expression_plot_grid_list = OUTPUT_PATH + "Rdata/LMM-DEG_volcano-plot_grids.rds",
        DE_genes_table = OUTPUT_PATH + "tabular/LMM-differential-expression_results.csv"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/differential_expression_analysis.R",
        output_path = OUTPUT_PATH,
        current_module = "differential_expression_analysis",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/differential-expression-analysis.out",
        err = OUTPUT_PATH + "logs/differential-expression-analysis.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} 1> {log.out} 2> {log.err}"

rule unsupervised_analysis:
    input:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_unsupervised-analysis.rds",
        R_file_unsupervised_clustering_plot_list = OUTPUT_PATH + "Rdata/plot-list_unsupervised-clustering.rds",
        R_file_unsupervised_clustering_plot_grid_list = OUTPUT_PATH + "Rdata/plot-list_unsupervised-clustering-grids.rds",
        R_file_16s_score_plot_list = OUTPUT_PATH + "Rdata/plot-list_16s-score.rds",
        R_file_unsupervised_clustering_cv_heatmap_list = OUTPUT_PATH + "Rdata/plot-list_cv_heatmaps.rds",
        CV_file = OUTPUT_PATH + "tabular/CV_results_by-normalization-method.csv"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/unsupervised_analysis.R",
        output_path = OUTPUT_PATH,
        current_module = "unsupervised_analysis",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/unsupervised-analysis.out",
        err = OUTPUT_PATH + "logs/unsupervised-analysis.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} 1> {log.out} 2> {log.err}"

rule analysis_16s: 
    input:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_normalized.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_16S-analysis.rds"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/16S_analysis.R",
        output_path = OUTPUT_PATH,
        current_module = "16S_analysis",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/16S-analysis.out",
        err = OUTPUT_PATH + "logs/16S-analysis.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} 1> {log.out} 2> {log.err}"

rule normalization:
    input:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-probes.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_normalized.rds",
        R_file_normalization_plot_list = OUTPUT_PATH + "Rdata/normalization_plot_list.rds"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/normalization.R",
        output_path = OUTPUT_PATH,
        current_module = "normalization",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/normalization.out",
        err = OUTPUT_PATH + "logs/normalization.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} 1> {log.out} 2> {log.err}"

rule qc_probes:
    input:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-segments.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-probes.rds",
        R_file_raw_qc_dat = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-probes-raw.rds",
        R_file_probe_qc_plot_list = OUTPUT_PATH + "Rdata/qc-probes_plot_list.rds",
        R_file_probe_qc_table = OUTPUT_PATH + "Rdata/qc-probes_table.rds"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/qc_probes.R",
        output_path = OUTPUT_PATH,
        current_module = "qc_probes",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/qc-probes.out",
        err = OUTPUT_PATH + "logs/qc-probes.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} 1> {log.out} 2> {log.err}"

rule qc_segments:
    input:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-study-design.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-segments.rds",
        R_file_main_module = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-segments_main-module.rds",
        R_file_segment_qc_summary_table = OUTPUT_PATH + "Rdata/qc-segments_summary_table.rds",
        R_file_segment_qc_plot_list = OUTPUT_PATH + "Rdata/qc-segments_plot_list.rds",
        Shiny_app = OUTPUT_PATH + "qc_probes_shiny_app.R"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/qc_segments.R",
        output_path = OUTPUT_PATH,
        current_module = "qc_segments",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/qc-segments.out",
        err = OUTPUT_PATH + "logs/qc-segments.err" 
    shell:
        """
        Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} 1> {log.out} 2> {log.err}
        cp src/qc_probes_shiny_app.R {params.output_path}
        cp src/qc_probes_protein_shiny_app.R {params.output_path}
        """

rule qc_study_design:
    input:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_raw.rds"
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_qc-study-design.rds",
        R_file_main_module = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_raw_main-module.rds",
        R_file_pkc_summary_table = OUTPUT_PATH + "Rdata/pkc_summary_table.rds",
        Shiny_app = OUTPUT_PATH + "qc_segments_shiny_app.R"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/qc_study-design.R",
        output_path = OUTPUT_PATH,
        current_module = "qc_study_design",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/qc-study-design.out",
        err = OUTPUT_PATH + "logs/qc-study-design.err" 
    shell:
        """
        Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} {input.R_file} 1> {log.out} 2> {log.err}
        cp src/qc_segments_shiny_app.R {params.output_path}
        """

rule data_import_cleaning:
    output:
        R_file = OUTPUT_PATH + "Rdata/NanoStringGeoMxSet_raw.rds"
    params:
        workflow_system = WORKFLOW_SYSTEM,
        script = "src/data_import_cleaning.R",
        output_path = OUTPUT_PATH,
        current_module = "data_import_cleaning",
        config_path = CONFIG_PATH
    log:
        out = OUTPUT_PATH + "logs/data-import-cleaning.out",
        err = OUTPUT_PATH + "logs/data-import-cleaning.err" 
    shell:
        "Rscript {params.script} {params.config_path} {params.workflow_system} {params.current_module} {params.output_path} 1> {log.out} 2> {log.err}"