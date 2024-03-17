
<h1 align="center">
  <br>
  <!--<a href="http://www.amitmerchant.com/electron-markdownify"><img src="https://raw.githubusercontent.com/amitmerchant1990/electron-markdownify/master/app/img/markdownify.png" alt="Markdownify" width="200"></a>-->
  <br>
  bisantine-geo
  <br>
</h1>

(<b>B</b>ugs <b>i</b>n <b>S</b>pace <b>An</b>alysis <b>T</b>oolkit <b>In</b>tegrated for <b>Efficiency</b>) — a Snakemake pipeline for processing NanoString GeoMx data.

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
  <a href="#credits">Credits</a> •
  <a href="#license">License</a>
</p>

## About

bisantine-geo is a Snakemake-powered pipeline that takes DCC files from a GeoMx experiment and performs quality control, normalization, basic unsupervised analysis, differential expression, and pathway analysis on your spatial transcriptomics data. 

## Usage
### Installing software requirements
To clone and run this pipeline, you'll need to have the following software installed on your machine:
1) [git](https://git-scm.com)
2) Some kind of conda package manager—e.g. [Anaconda](https://www.anaconda.com/download), [miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/), or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). I personally use miniconda.
3) Snakemake. Once you have your conda package manager set up, I (and the Snakemake people) recommend [using that to install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Cloning the repository
A.k.a. downloading the pipeline to your machine. 
On the command line, navigate to the directory in which you want to download the pipeline, then use `git clone` to download it.
```bash
# Navigate to the directory where you want to save the pipeline.
cd path/to/directory/

# Clone the repository.
$ git clone https://github.com/alonzowolfram/bisantine-geo
```
> **Note**
> Replace `path/to/directory/` with the path to the target directory on your machine.

### Configuration
Now that we have the pipeline downloaded, we will need to configure the settings for your particular experiment. 

All the settings you will need to edit are stored in one convenient file, `config.yaml`, located in the `bisantine-geo` folder you just downloaded. Open this file in your favorite text editor and edit the settings accordingly. See the comments above each setting for documentation.

### Running the pipeline


## Credits

bisantine-geo is developed and maintained by the Digital Spatial Profiling team of the PRIME-TR platform at the University of Texas MD Anderson Cancer Center.

## License

MIT
