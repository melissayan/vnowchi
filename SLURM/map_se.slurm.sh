#!/usr/bin/bash
#
#SBATCH --job-name=map_se
#SBATCH --ntasks=1
#SBATCH --get-user-env
#SBATCH --cpus-per-task=12
#SBATCH --mem=48000
#SBATCH --partition=highio
#SBATCH --time=2160
#SBATCH --output=/results/SLURM_debug/map_se_%A_%a.out
#SBATCH --error=/results/SLURM_debug/map_se_%A_%a.out

set -e
set -x


#----------------------------------------------------------------------------------------
# Map with BWA_MEM for single-end reads
#----------------------------------------------------------------------------------------
SCRIPTS_DIR=$1
FILES=$2
PATH=$3
WORKING_DIR=$4
REF=$5

FILELINE=$((SLURM_ARRAY_TASK_ID + 1 ))
FILE=$(sed "${FILELINE}q;d" "$FILES")

srun $SCRIPTS_DIR/run_bwa_mem_se.sh $FILE $PATH $WORKING_DIR $REF


