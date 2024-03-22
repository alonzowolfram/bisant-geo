
<h1 align="center">
  <br>
  <!--<a href="http://www.amitmerchant.com/electron-markdownify"><img src="https://raw.githubusercontent.com/amitmerchant1990/electron-markdownify/master/app/img/markdownify.png" alt="Markdownify" width="200"></a>-->
  <br>
  bisantine-geo
  <br>
</h1>

(<b>B</b>ugs <b>i</b>n <b>S</b>pace <b>An</b>alysis <b>T</b>oolkit <b>In</b>tegrated for <b>E</b>fficiency) — a Snakemake pipeline for processing NanoString GeoMx data.

<!--
<p align="center">
  <a href="https://badge.fury.io/js/electron-markdownify">
    <img src="https://badge.fury.io/js/electron-markdownify.svg"
         alt="Gitter">
  </a>
  <a href="https://gitter.im/amitmerchant1990/electron-markdownify"><img src="https://badges.gitter.im/amitmerchant1990/electron-markdownify.svg"></a>
  <a href="https://saythanks.io/to/bullredeyes@gmail.com">
      <img src="https://img.shields.io/badge/SayThanks.io-%E2%98%BC-1EAEDB.svg">
  </a>
  <a href="https://www.paypal.me/AmitMerchant">
    <img src="https://img.shields.io/badge/$-donate-ff69b4.svg?maxAge=2592000&amp;style=flat">
  </a>
</p>
-->

<p align="center">
  <a href="#about">About</a> •
  <a href="#usage">Usage</a> •
  <a href="#changelog">Changelog</a> •
  <a href="#credits">Credits</a> •
  <a href="#license">License</a>
</p>

## About

bisantine-geo is a Snakemake-powered pipeline that takes DCC files from a GeoMx experiment and performs quality control, normalization, basic unsupervised analysis, differential expression, and pathway analysis on your spatial transcriptomics data. 

## Usage
### Installing software requirements
To clone and run this pipeline, you'll need to have the following software installed on your machine:
1) [git](https://git-scm.com)
2) Some kind of conda package manager—e.g. [Anaconda](https://www.anaconda.com/download), [miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/), or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). I personally use mamba, for its speed.
3) Snakemake. Once you have your conda package manager set up, I (and the Snakemake people) recommend [using that to install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Cloning the repository
A.k.a. downloading the pipeline to your machine. 
On the command line, navigate to the directory in which you want to download the pipeline, then use `git clone` to download it.
```bash
# Navigate to the directory where you want to save the pipeline.
cd path/to/directory/

# Clone the repository.
git clone https://github.com/alonzowolfram/bisantine-geo
```
> **Note:**
> Replace `path/to/directory/` with the path to the target directory on your machine.

### Configuration
Now that you have the pipeline downloaded, you will need to configure the settings for your particular experiment. 

All the settings you will need to edit are stored in one convenient file, `config.yaml`, located in the `bisantine-geo` folder you just downloaded. Open this file in your favorite text editor and edit the settings accordingly. See the comments above each setting for documentation. You can even store different configuration files in the `profiles` folder (or anywhere, really—it's just there for convenience) and call them with the `--configfiles` flag to the `snakemake` command when running the pipeline (see the section "Running the pipeline," below.)

### Setting up your conda environment
Using the `environment.yaml` file in the `bisantine-geo/envs` folder, create a conda environment. See directions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) on how to set up a conda environment from a YAML file. 

> **Note:**
> If you're on a Mac with an x86_64 ISA, you might have to prefix your conda environment creation command with `CONDA_SUBDIR=osx-64`. <br />
> For example, `CONDA_SUBDIR=osx-64 mamba env create -n geomx -f environment.yml`

### Running the pipeline
Now that you've set up your conda environment and the configuration file for your run, you can run the pipeline. To run bisantine-geo, activate the conda environment you created above, then navigate into the `bisantine-geo` folder. From there, running Snakemake is quite simple:

```bash
NCORES=3
snakemake --configfiles path/to/config.yaml --cores $NCORES
```

Replace `3` with the desired number of cores. If you're familiar with Snakemake, you can run individual modules as well, saving you time if you want to pick up from a previously completed run. For example:

```bash
snakemake --configfiles path/to/config.yaml --cores $NCORES -R --until unsupervised_analysis
```
will run the pipeline up to the "unsupervised_analysis" module.
<br />
<br />
A tutorial with more details is forthcoming. Sometime. Maybe. 
<br />
<br />
After the pipeline is finished running, the folder containing the output files will be available in the `out` directory. See the `README.md` file in the `out` directory for an explanation of the outputs. 

## Changelog
<b>2024/03/22</b> - v0.2.0-alpha:
* Lots and lots of bug fixes.
* Added output: list of PKC modules.
* Added output: NTC summary.
* Changed default values of config.YAML to be more indicative of what should go there.
* Added ability to specify YAML config file via command line; added `profiles` folder for storage of YAML config files.
* Added Shiny app to see how different parameter changes affect number of viable segments.
* Cutoffs on segment QC graphs are now based on user-input values.
* Added functionality to correct names in annotation sheet to match functions' expected input.
* Added ability to have multiple compartment variables, including variables that are combinations of two or more existing variables.
* Added ability to start pipeline from a previous run. 
* Added ability for user to manually enter the PKC files to be used.

<b>2024/03/17</b> - v0.1.0-alpha released on GitHub.

## Credits

bisantine-geo is developed and maintained by the Digital Spatial Profiling team of the PRIME-TR platform at the University of Texas MD Anderson Cancer Center. Some of the code is adapted from the GeomxTools Bioconductor package vignette, found [here](https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html).

## License

bisantine-geo is released under the GNU General Public License v3.0. For more information, see the `LICENSE` file. 
