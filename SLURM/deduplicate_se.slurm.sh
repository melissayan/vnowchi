#!/usr/bin/bash
#
#SBATCH --job-name=dedup_se
#SBATCH --ntasks=1
#SBATCH --get-user-env
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --partition=highio
#SBATCH --time=2160
#SBATCH --output=/results/SLURM_debug/dedup_se_%A_%a.out
#SBATCH --error=/results/SLURM_debug/dedup_se_%A_%a.out

set -e
set -x


#----------------------------------------------------------------------------------------
# Trim low quality bases w/ Trimmomatic for single-end reads
#----------------------------------------------------------------------------------------
SCRIPTS_DIR=$1
FILES=$2
PATH=$3

FILELINE=$((SLURM_ARRAY_TASK_ID + 1 ))
FILE=$(sed "${FILELINE}q;d" "$FILES")

srun $SCRIPTS_DIR/run_fastx_collapser_se.sh $FILE $PATH


