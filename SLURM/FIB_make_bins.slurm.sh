#!/usr/bin/bash
#
#SBATCH --job-name=make_bins
#SBATCH --ntasks=1
#SBATCH --get-user-env
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --partition=highio
#SBATCH --time=2160
#SBATCH --output=/results/SLURM_debug/make_bins_%A.out
#SBATCH --error=/results/SLURM_debug/make_bins_%A.out

set -e
set -x


#----------------------------------------------------------------------------------------
# Make bins from fibroblasts
#----------------------------------------------------------------------------------------
SCRIPTS_DIR=$1
BAM_FILE=$2
COUNT=$3
OUTPUT=$4

srun $SCRIPTS_DIR/make_bins.sh $SCRIPTS_DIR $BAM_FILE $COUNT $OUTPUT


