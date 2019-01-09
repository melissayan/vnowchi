#!/usr/bin/bash
#
#SBATCH --job-name=get_cnv
#SBATCH --ntasks=1
#SBATCH --get-user-env
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --partition=highio
#SBATCH --time=2160
#SBATCH --output=/results/SLURM_debug/get_cnv_%A_%a.out
#SBATCH --error=/results/SLURM_debug/get_cnv_%A_%a.out

set -e
set -x


#----------------------------------------------------------------------------------------
# Get copy number and classify samples
#----------------------------------------------------------------------------------------
SCRIPTS_DIR=$1
FILES=$2
CHR_SIZE_TABLE=$3
SUMMARY_FILE=$4

FILELINE=$((SLURM_ARRAY_TASK_ID + 1 ))
FILE=$(sed "${FILELINE}q;d" "$FILES")

srun $SCRIPTS_DIR/get_copy_number_and_classify.R $FILE $CHR_SIZE_TABLE $SUMMARY_FILE


