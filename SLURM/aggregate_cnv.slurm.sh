#!/usr/bin/bash
#
#SBATCH --job-name=aggregate
#SBATCH --ntasks=1
#SBATCH --get-user-env
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --partition=highio
#SBATCH --time=2160
#SBATCH --output=/results/SLURM_debug/aggregate_%A_%a.out
#SBATCH --error=/results/SLURM_debug/aggregate_%A_%a.out

set -e
set -x


#----------------------------------------------------------------------------------------
# Aggregate sample data by embryo
#----------------------------------------------------------------------------------------
SCRIPTS_DIR=$1
FILES=$2
GROUP_INFO_TABLE=$3

FILELINE=$((SLURM_ARRAY_TASK_ID + 1 ))
FILE=$(sed "${FILELINE}q;d" "$FILES")

srun $SCRIPTS_DIR/aggregate_cnv.R $FILE $GROUP_INFO_TABLE 


