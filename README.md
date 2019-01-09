# Variable Non-Overlapping Window CBS and HMM Intersect (VNOWCHI): A Copy Number Variant Calling Pipeline
-----------------------------------------------------------------------------------------
[![R: v3.5.0](https://img.shields.io/badge/R-v3.5.0-198ce7.svg)](https://cran.rstudio.com/)  [![DOI:10.1101/241851 ](https://img.shields.io/badge/DOI-10.110/241851-70db70.svg)](https://doi.org/10.1101/241851)

VNOWCHI is a copy number variant calling pipeline that utilizes the Circular Binary Segmentaion and Hidden Markov Model algorithms to determine the ploidy and sex status of embryos.  

Detailed information regarding VNOWCHI and the CBS and HMM algorithms can be found on the [Wiki](https://github.com/melissayan/vnowchi/wiki/Algorithm-Information).

## Table of Contents
* [Getting Started](#getting-started)
	* [Prerequisites](#prerequisites)
	* [Installation](#installation)
	* [How to Use](#how-to-use)
* [Results](#results)
* [Additonal Notes](#additional-notes)
* [Authors](#authors)
* [Acknowledgments](#acknowledgements)
* [Citation](#citation)
* [License](#license)

## Getting Started
The scripts provided are tailored for use on [OHSU's exacloud](https://www.ohsu.edu/xd/research/research-cores/advanced-computing-center/exacloud.cfm) server using SLURM.  For help with SLURM, please refer to ACC's [tutorials](https://accdoc.ohsu.edu/exacloud/guide/getting-started/).

To use scripts without SLURM, comment `srun` and `sbatch` commands and uncomment the commands beneath.  

### Prerequisites
* Linux environment
* The following software tools installed:
    * FastQC 0.10.1 
    * Trimmomatic 0.35
    * FASTX-Toolkit 0.0.13
    * BWA-MEM 0.7.9a-r786
    * SAMtools 0.1.19-44428cd
    * BEDtools 2.25.0
    * FastUniq 1.1
* R version 3.5.0 with the following R packages installed:
    * GenomicAlignments
    * DNAcopy
    * HMMcopy
    * IRanges
    * GenomicRanges
    * dplyr
    * ggplot2
