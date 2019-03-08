![VNOWCHI](https://github.com/melissayan/vnowchi/blob/master/img/vnowchi.png)
-----------------------------------------------------------------------------------------
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


## Getting Started
The scripts provided are tailored for use on [OHSU's exacloud](https://www.ohsu.edu/xd/research/research-cores/advanced-computing-center/exacloud.cfm) server using SLURM.  For help with SLURM, please refer to ACC's [tutorials](https://accdoc.ohsu.edu/exacloud/guide/getting-started/).

To use scripts without SLURM, comment `srun` and `sbatch` commands and uncomment the commands beneath.  

### Prerequisites
* Linux environment
* Java 8 ([How to install Java](https://github.com/in28minutes/getting-started-in-5-steps)) for Trimmomatic
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

### Installation
Installation instructions can be found on the [Wiki Installation page](https://github.com/melissayan/vnowchi/wiki/Installation).

### How to Use
Detailed "How to Use" instructions are located on the [Wiki How to Use page](https://github.com/melissayan/vnowchi/wiki/How-to-Use).
#### 1. Modify scripts with your directory information.
#### 2. Copy fibroblast FASTQs and sample FASTQs to the following locations:
```
Fibroblast samples (5 scDNA-seq samples preferred):
/your/working/dir/CopyNumberPipeline/results/FIBROBLASTS/FASTQ

Single-ended samples:
/your/working/dir/CopyNumberPipeline/results/fastq/SE

Paired-ended samples:
/your/working/dir/CopyNumberPipeline/results/fastq/PE
```
#### 3. Generate bins for the pipeline using FIBROBLAST data.
```
sbatch PIPELINE_bins.sh
```
#### 4. Run pipeline on actual data.
```
sbatch PIPELINE_VNOWCHI.sh
```
#### 5. Look at results:
* CNV plots for all samples, includes all and individual chromosomes
* Mapping summary statistics for VNOWCHI_summary.txt
* Tabluar summary for all samples classified by ploidy and sex status
* Tabular summary for all embryos classified by ploidy and sex status based on samples

## Results:
* Mapping summary statistics `VNOWCHI_Summary.txt`
* Tabular summary of all individual samples CNV calls with embryo, blastomere, ploidy and sex classifications
	* `CNV_<SE|PE>_<bin>.sampleSummary.txt`
* Tabular summary of all embryos classified by ploidy and sex status
	* `CNV_<SE|PE>_<bin>.embryoSummary.txt`
* CNV plots for all samples by chromosome or by all chromosomes 
	* `<sampleName>_<chromosome>.png`
	* `<sampleName>_<all>.png`
	
<kbd><img src="https://github.com/melissayan/vnowchi/blob/master/img/samplePlot.png"></kbd>

## Additional Notes 
* Please note that R package dply will behave differently than intended if R package plyr is loaded. More info regarding the issue can be found [here](https://github.com/tidyverse/dplyr/issues/29) and [here](https://github.com/tidyverse/dplyr/issues/347).  Here's a possible solution from [Stack Overflow](https://stackoverflow.com/questions/22801153/dplyr-error-in-n-function-should-not-be-called-directly) if you get any errors. 
* Might need to modify step 3 in `PIPELINE_bins.sh` and `PIPELINE_VNOWCHI.sh` to ensure script will accept the provided fastq file name format pattern
* If trying to use Rscripts in Rstudio, some scripts have issues.  Ex. get_copy_number.R has no problems on server but does not work in RStudio, could be related to R version. 

## Authors
* **Melissa Yan** - extended the VNOWC pipeline to include CHI, classify samples/embryos, accommodate different genomes, and run on SLURM
* [**Nathan Lazar**](https://github.com/nathanlazar) - original author of [Variable Non-Overlapping Window CBS (VNOWC)](https://github.com/nathanlazar/Oocyte_CN)
* **Kristof Torkency** - original author of CBS/HMM Intersect (CHI) pipeline

## Acknowledgments
This project would not be possible without the support from the Chavez Lab, Carbone Lab, Adey Lab, and the Biostatistics & Bioinformatics Core:
* [Chavez Lab](https://www.ohsu.edu/people/shawnchavez/afe058a49a4aac1f07559cb5b32c424a):
	* Brittany L. Daughtry
	* Kelsey E. Brooks
	* Jimi L. Rosenkrantz
	* Shawn L. Chavez
* [Carbone Lab](http://carbonelab.com/): 
	* [Nathan H. Lazar](https://github.com/nathanlazar)
	* Brett Davis
	* Lucia Carbone
* [Adey Lab](https://adeylab.org/andrew-adey/):
	* Kristof A. Torkenczy
	* Andrew Adey
* [Biostatistics & Bioinformatics Core](http://www.ohsu.edu/bbc)
	* Suzi S. Fei

## Citation
Daughtry, B. L., Rosenkrantz, J. L., Lazar, N. H., Fei, S. S., Redmayne, N., Torkenczy, K. A., Adey, A., Gao, L., Park, B., Nevonen, K.A., Carbone, L., Chavez, S. L. (2019). Single-cell sequencing of primate preimplantation embryos reveals chromosome elimination via cellular fragmentation and blastomere exclusion. *Genome Research*, *29(3)*, 367-382. [doi:10.1101/gr.239830.118](https://genome.cshlp.org/content/early/2019/01/25/gr.239830.118)
