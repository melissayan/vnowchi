#!/usr/bin/bash
#
#SBATCH --job-name=FIB_trim_se
#SBATCH --ntasks=1
#SBATCH --get-user-env
#SBATCH --cpus-per-task=12
#SBATCH --mem=4000
#SBATCH --partition=highio
#SBATCH --time=2160
#SBATCH --output=/results/SLURM_debug/FIB_trim_se_%A_%a.out
#SBATCH --error=/results/SLURM_debug/FIB_trim_se_%A_%a.out

set -e
set -x


#----------------------------------------------------------------------------------------
# Trim low quality bases w/ Trimmomatic for single-end FIBROBLAST reads
#----------------------------------------------------------------------------------------
trim_dir=$1
OUT_DIR=$2
FILES=$3

FILELINE=$((SLURM_ARRAY_TASK_ID + 1 ))
F=$(sed "${FILELINE}q;d" "$FILES")

# Might need to change this based on the fastq filename
name=${F/_001.fastq.gz/}; name=${F/.gz/}; name=${name/.fastq.gz}; name=$(basename $name)


java -jar $trim_dir/trimmomatic-0.35.jar \
    SE -threads 12 -phred33 $F \
    ${OUT_DIR}/${name}.trim.fq.gz  \
    ILLUMINACLIP:$trim_dir/adapters/TruSeq3-SE.fa:2:30:10 LEADING:5 TRAILING:5 \
    SLIDINGWINDOW:4:15 MINLEN:36 &> TRIM_SE/${name}_trim.txt

