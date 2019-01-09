#!/usr/bin/bash
#
#SBATCH --get-user-env
#SBATCH --partition=exacloud
#SBATCH --cpus-per-task=12
#SBATCH --mem=48000
#SBATCH --time=2160
#SBATCH --verbose
#SBATCH --output=/results/SLURM_debug/bins_%j.out
#SBATCH --error=/results/SLURM_debug/bins_%j.out
#SBATCH --export=ALL  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

set -e
set -x


#----------------------------------------------------------------------------------------
# Constructs Bins Using Fibroblast Samples
#----------------------------------------------------------------------------------------
# Steps to make variable sized bins using reads from fibroblasts known to be 
# euploid. We only use the single-cell samples to construct bins
# Bins are made treating reads as paired-end and single-ended with 
# 500, 1,000, 2,000 and 4,000 reads per bin
#
# Software versions:
#   - FastQC: v0.10.1
#   - Trimmomatic v0.35  *note: .jar file needs to be in specific directory see trim_dir
#   - FASTX-Toolkit: v0.0.13
#   - BWA-MEM: 0.7.9a-r786
#   - SAMtools: 0.1.19-44428cd
#   - BEDtools: v2.25.0
#   - FastUniq: v1.1
#
# Modified by Melissa Yan for SLURM, to include concatenating single cell XX paired-end
# read sets and summary stats.
# 
# Built on original script by Nathan Lazar.
#   - Original script:
#       Project Title: Oocyte_CN
#       Script Title: PIPELINE_bins
#       Author: Nathan Lazar (nathan dot lazar at gmail dot com)
#       Date: September 19, 2017
#       Availability: https://github.com/nathanlazar/Oocyte_CN/blob/master/PIPELINE_bins
#----------------------------------------------------------------------------------------


CURRENT_DIR=$(pwd)
RESULTS_DIR=$(pwd)/results/
WORKING_DIR=$(pwd)/results/FIBROBLASTS
SCRIPTS_DIR=$(pwd)/scripts
SLURM_DIR=$(pwd)/SLURM
SLURM_DEBUG=$WORKING_DIR/SLURM_DEBUG
TEMP_DIR=$(pwd)/temp

#change this to your username
USERNAME=

# Add paths to tools
PATH=$PATH/:/home/users/$USERNAME/bin
PATH=$PATH/:/home/users/$USERNAME/bin/bedtools2/bin
PATH=$PATH/:/home/users/$USERNAME/bin/FastQC
PATH=$PATH/:/home/users/$USERNAME/bin/fastx/bin
PATH=$PATH/:/home/users/$USERNAME/bin/FastUniq/source
export PATH
trim_dir=/home/users/$USERNAME/bin/Trimmomatic-0.35

# Set the reference genome location
REF=

# Set the number of parallel threads to use
cores=12

#----------------------------------------------------------------------------------------
# 0) Collect all raw read files to process by copying fibroblasts to fastq directory
# *note: only use single cell XX paired-end readsets
#
# Controls: xx paired-end read sets,
# 1 with 0 cells:			170725_Chavez.BOV_FIB_NTC
# 5 with 1 cell:			170725_Chavez.BOV_FIB_XX_1a - 170725_Chavez.BOV_FIB_XX_1e
# 1 with 5 cells:			170725_Chavez.BOV_FIB_XX_5
# 1 with 10 cells:			170725_Chavez.BOV_FIB_XX_10
# 1 with all XX paired-end combined**:	xx-euploid_all
# 1 with just single cells combined**:	xx-euploid-1_all
#
# **Generate last two just by concatentating the unzipped files in step 1 below
#----------------------------------------------------------------------------------------
if [ ! -d $TEMP_DIR ]; then mkdir $TEMP_DIR; fi
cd $WORKING_DIR

<<"COMMENT"
#COMMENT
#----------------------------------------------------------------------------------------
# 1. Concatenate single cell XX paired-end read sets 
# *note: May have to change this part based on fastq filename patterns
#       *xx_1[a-z]*     *xx_1[letter from a-z]
#       *XX_1[a-z]*     
#----------------------------------------------------------------------------------------
cd FASTQ
touch xx-euploid_all.R1; touch xx-euploid_all.R2
touch xx-euploid-1_all.R1; touch xx-euploid-1_all.R2
for f in `find . \( -name "*xx*" -o -name "*XX*" \)`; do
  if [[ $f != *"xx-euploid"* ]]; then  
    if [[ $f == *"R1"* ]]; then
      zcat $f >> xx-euploid_all.R1
    elif [[ $f == *"R2"* ]]; then
      zcat $f >> xx-euploid_all.R2
    fi
  fi
done
for f in `find . \( -name "*xx_1[a-z]*" -o -name "*XX_1[a-z]*" \)`; do
  if [[ $f != *"xx-euploid"* ]]; then      
    if [[ $f == *"R1"* ]]; then
      zcat $f >> xx-euploid-1_all.R1
    elif [[ $f == *"R2"* ]]; then
      zcat $f >> xx-euploid-1_all.R2
    fi
  fi
done
gzip xx-euploid*
cd ..

echo ------------------- done concatenate single cells;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 2. Run fastqc on the raw reads
#----------------------------------------------------------------------------------------
if [ ! -d FASTQC ]; then mkdir FASTQC; fi

FILES=$TEMP_DIR/FASTQ.txt
ls -1 $WORKING_DIR/FASTQ/*gz > $FILES
NUMFILES=$(($(cat $FILES | wc -l) - 1))
if [ $NUMFILES -ge 0 ]; then
  sbatch --array=0-$NUMFILES $SLURM_DIR/fastqc.slurm.sh $WORKING_DIR/FASTQC $FILES
fi

exit
#COMMENT
#***** Please check jobs completed before continuing.

# Look at the number of reads filtered, read length, etc. before trimming and write to 
# summary.txt
cd FASTQC
for f in `ls *.zip`; do unzip $f; done
grep "Total Sequences" */fastqc_data.txt > summary.txt
grep "Sequence length" */fastqc_data.txt >> summary.txt
grep "Total Deduplicated Percentage" */fastqc_data.txt >> summary.txt
grep "^%GC" */fastqc_data.txt >> summary.txt
rm *.html
rm -r *fastqc
cd ..

echo ------------------- done FASTQC on raw reads;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 3. Trim low quality bases w/ Trimmomatic
#----------------------------------------------------------------------------------------
if [ ! -d TRIM_SE ]; then mkdir TRIM_SE; fi
if [ ! -d TRIM_PE ]; then mkdir TRIM_PE; fi

# Trim reads as if they were single-ended
FILES=$TEMP_DIR/FASTQ_SE.txt
ls -1 $WORKING_DIR/FASTQ/*gz > $FILES
NUMFILES=$(($(cat $FILES | wc -l) - 1))
if [ $NUMFILES -ge 0 ]; then
  sbatch --array=0-$NUMFILES $SLURM_DIR/FIB_trim_se.slurm.sh $trim_dir $WORKING_DIR/TRIM_SE $FILES
fi

# Trim paired end reads
FILES=$TEMP_DIR/FASTQ_PE.txt
ls -1 $WORKING_DIR/FASTQ/*R1*gz > $FILES
NUMFILES=$(($(cat $FILES | wc -l) - 1))
if [ $NUMFILES -ge 0 ]; then
  sbatch --array=0-$NUMFILES $SLURM_DIR/FIB_trim_pe.slurm.sh $trim_dir $WORKING_DIR/TRIM_PE $FILES
fi

echo ------------------- done TRIM_PE and TRIM_SE;echo `date +"%Y-%m-%d %T"`;echo;

exit
#COMMENT
#***** Please check jobs completed before continuing.


#----------------------------------------------------------------------------------------
# 4. Run fastqc on the trimmed reads
#----------------------------------------------------------------------------------------
if [ ! -d FASTQC_TRIM_SE ]; then mkdir FASTQC_TRIM_SE; fi
if [ ! -d FASTQC_TRIM_PE ]; then mkdir FASTQC_TRIM_PE; fi

FILES=$TEMP_DIR/TRIM_SE.txt
ls -1 $WORKING_DIR/TRIM_SE/*trim.fq.gz > $FILES
NUMFILES=$(($(cat $FILES | wc -l) - 1))
if [ $NUMFILES -ge 0 ]; then
  sbatch --array=0-$NUMFILES $SLURM_DIR/fastqc.slurm.sh $WORKING_DIR/FASTQC_TRIM_SE $FILES
fi

FILES=$TEMP_DIR/TRIM_PE.txt
ls -1 $WORKING_DIR/TRIM_PE/*trim.fq.gz > $FILES
NUMFILES=$(($(cat $FILES | wc -l) - 1))
if [ $NUMFILES -ge 0 ]; then
  sbatch --array=0-$NUMFILES $SLURM_DIR/fastqc.slurm.sh $WORKING_DIR/FASTQC_TRIM_PE $FILES
fi

exit
#COMMENT
#***** Please check jobs completed before continuing.

# Look at the number of reads filtered, read length, etc. before trimming
# Write to summary.txt
cd FASTQC_TRIM_SE
for f in `ls *.zip`; do unzip $f; done
grep "Total Sequences" */fastqc_data.txt > summary.txt
grep "Sequence length" */fastqc_data.txt >> summary.txt
grep "Total Deduplicated Percentage" */fastqc_data.txt >> summary.txt
grep "^%GC" */fastqc_data.txt >> summary.txt
rm *.html; rm -r *fastqc; cd ..
cd FASTQC_TRIM_PE
for f in `ls *.zip`; do unzip $f; done
grep "Total Sequences" */fastqc_data.txt > summary.txt
grep "Sequence length" */fastqc_data.txt >> summary.txt
grep "Total Deduplicated Percentage" */fastqc_data.txt >> summary.txt
grep "^%GC" */fastqc_data.txt >> summary.txt
rm *.html; rm -r *fastqc; cd ..

echo ------------------- done FASTQC_TRIM_PE FASTQC_TRIM_SE;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 5. Deduplicate (a.k.a. collapse) reads that are exact copies before mapping
#    with fastx for single-end reads and fastuniq for paired-end reads
#----------------------------------------------------------------------------------------
if [ ! -d COLLAPSED_SE ]; then mkdir COLLAPSED_SE; fi
if [ ! -d COLLAPSED_PE ]; then mkdir COLLAPSED_PE; fi

# Collapse single-end reads
FILES=$TEMP_DIR/trim_SE.txt
ls -1 $WORKING_DIR/TRIM_SE/*.trim.fq.gz > $FILES
NUMFILES=$(($(cat $FILES | wc -l) - 1))
if [ $NUMFILES -ge 0 ]; then
  sbatch --array=0-$NUMFILES $SLURM_DIR/deduplicate_se.slurm.sh $SCRIPTS_DIR $FILES $PATH
fi

# Collapse paired-end reads
FILES=$TEMP_DIR/trim_PE.txt
ls -1 $WORKING_DIR/TRIM_PE/*R1*.trim.fq.gz > $FILES
NUMFILES=$(($(cat $FILES | wc -l) - 1))
if [ $NUMFILES -ge 0 ]; then
  sbatch --array=0-$NUMFILES $SLURM_DIR/deduplicate_pe.slurm.sh $SCRIPTS_DIR $FILES $PATH
fi

exit
#COMMENT
#***** Please check jobs completed before continuing.

echo ------------------- done COLLAPSED_SE COLLAPSED_PE;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 6. Map with BWA-MEM  
#    *note: The genome must first be indexed w/ bwa index. Example for rheMac8 genome:
#          REF=$WORKING_DIR/RheMac8/rheMac8.fa
#          name=$(basename $REF); name=$WORKING_DIR/RheMac8/$name;
#          bwa index -a bwtsw $name $name.amb $name.ann $name.bwt $name.pac $name.sa
#
#    This script runs bwa-mem, and filters for quality of 30
#    Outputs are:
#      <read_name>.mapped_reads   # Read counts at each step
#      <read_name>.q30.sort.bam
#      <read_name>.q30.sort.bam.bai
#----------------------------------------------------------------------------------------
if [ ! -d MAPPED_SE ]; then mkdir MAPPED_SE; fi
if [ ! -d MAPPED_PE ]; then mkdir MAPPED_PE; fi

# Map single-end reads
FILES=$TEMP_DIR/map_SE.txt
ls -1 $WORKING_DIR/COLLAPSED_SE/*.fq.gz > $FILES
NUMFILES=$(($(cat $FILES | wc -l) - 1))
if [ $NUMFILES -ge 0 ]; then
  sbatch --array=0-$NUMFILES $SLURM_DIR/map_se.slurm.sh $SCRIPTS_DIR $FILES $PATH $WORKING_DIR $REF
fi

# Map paired-end reads
FILES=$TEMP_DIR/map_PE.txt
ls -1 $WORKING_DIR/COLLAPSED_PE/*.fq.gz > $FILES
NUMFILES=$(($(cat $FILES | wc -l) - 1))
if [ $NUMFILES -ge 0 ]; then
  sbatch --array=0-$NUMFILES $SLURM_DIR/map_pe.slurm.sh $SCRIPTS_DIR $FILES $PATH $WORKING_DIR $REF
fi

exit
COMMENT
#***** Please check jobs completed before continuing.

## TODO: there was an error mapping individual samples with paired ends

# Write mapped results for each sample after alignment to summary_MAPPED.txt:
#    - Mapped Reads
#    - Reads Passing Quality of 30
#    - Unique Mapping Positions
srun $SCRIPTS_DIR/summarize_v2.sh $WORKING_DIR
#$SCRIPTS_DIR/summarize_v2.sh $WORKING_DIR

echo ------------------- done MAPPED_SE MAPPED_PE and stats;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 7. Combine reads that were mapped as single ended and look at the read depths
#----------------------------------------------------------------------------------------
# Combine single-end reads
for f1 in `ls $WORKING_DIR/MAPPED_SE/*R1*.bam`; do
  name=${f1/.R1.q30.sort.bam}; f2=${f1/R1/R2}
  srun samtools merge ${name}_merge.bam $f1 $f2 &
done
wait

# Make a .txt file of the read depth at each position
for f in `ls $WORKING_DIR/MAPPED_SE/*.bam`; do
  srun samtools depth $f > ${f/.bam/.depth} &
done
for f in `ls $WORKING_DIR/MAPPED_PE/*.bam`; do
  srun samtools depth $f > ${f/.bam/.depth} &
done
wait

echo ------------------- done depth files;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 8. Convert bam files to bed files for viewing
#----------------------------------------------------------------------------------------
for f in `ls $WORKING_DIR/MAPPED_SE/*.bam`; do
  srun bamToBed -i $f > ${f/.bam/.bed} &
done
for f in `ls $WORKING_DIR/MAPPED_PE/*.bam`; do
  srun bamToBed -i $f > ${f/.bam/.bed} &
done
wait

echo ------------------- done bed files;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 9. Summarize mapping statistics w/ R
#----------------------------------------------------------------------------------------
# TODO: rework this script
srun ${SCRIPTS_DIR}/get_stats.R "$WORKING_DIR" FIBROBLASTS
#./${SCRIPTS_DIR}/get_stats.R "$WORKING_DIR" FIBROBLASTS

echo ------------------- done get_stats files;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 10. Make bins using the combined single-cell data with expected number of reads
#     of 500, 1,000, 2,000 and 4,000
#----------------------------------------------------------------------------------------
if [ ! -d BINS ]; then mkdir BINS; fi
if [ ! -d BINS/SE ]; then mkdir BINS/SE; fi
if [ ! -d BINS/PE ]; then mkdir BINS/PE; fi

for cnt in 500 1000 2000 4000; do
  sbatch ${SLURM_DIR}/FIB_make_bins.slurm.sh $SCRIPTS_DIR $WORKING_DIR/MAPPED_SE/xx-euploid-1_all_merge.bam $cnt $WORKING_DIR/BINS/SE/xx-euploid-1_all_merge_${cnt}.bins 
done
for cnt in 500 1000 2000 4000; do
  sbatch ${SLURM_DIR}/FIB_make_bins.slurm.sh $SCRIPTS_DIR $WORKING_DIR/MAPPED_PE/xx-euploid-1_all.q30.sort.bam $cnt $WORKING_DIR/BINS/PE/xx-euploid-1_all_${cnt}.bins
done
wait

echo ------------------- done make_bins stats;echo `date +"%Y-%m-%d %T"`;echo;


