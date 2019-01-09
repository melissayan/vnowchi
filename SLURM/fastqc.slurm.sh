#!/usr/bin/bash
#
#SBATCH --job-name=fastq
#SBATCH --ntasks=1
#SBATCH --get-user-env
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --partition=highio
#SBATCH --time=2160
#SBATCH --output=/results/SLURM_debug/fastq_%A_%a.out
#SBATCH --error=/results/SLURM_debug/fastq_%A_%a.out

set -e
set -x


#----------------------------------------------------------------------------------------
# Run FastQC on raw reads
#----------------------------------------------------------------------------------------
OUT_DIR=$1
FILES=$2

FILELINE=$((SLURM_ARRAY_TASK_ID + 1 ))
FILENAME=$(sed "${FILELINE}q;d" "$FILES")

fastqc --outdir $OUT_DIR --threads 12 $FILENAME


