#!/usr/bin/bash
#
#SBATCH --job-name=counts
#SBATCH --ntasks=1
#SBATCH --get-user-env
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --partition=highio
#SBATCH --time=2160
#SBATCH --output=/results/SLURM_debug/counts_%A_%a.out
#SBATCH --error=/results/SLURM_debug/counts_%A_%a.out

set -e
set -x


#----------------------------------------------------------------------------------------
# Get counts and ratios in the specified bin 
#----------------------------------------------------------------------------------------
SCRIPTS_DIR=$1
BIN_FILE=$2
FILES=$3
WORKING_DIR=$4

FILELINE=$((SLURM_ARRAY_TASK_ID + 1 ))
FILE=$(sed "${FILELINE}q;d" "$FILES")
out=$(basename $FILE); out=${out/.q30.sort.bam/.counts}

srun $SCRIPTS_DIR/count_bin_reads.R $BIN_FILE $FILE $WORKING_DIR/$out 


