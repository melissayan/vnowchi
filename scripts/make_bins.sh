#!/bin/bash

# Usage: make_bins.sh <reads.bam> <reads per bin> <out.bins>

#dir=`pwd`
#dir=${dir/FIBROBLASTS/}
#Rscript $dir/make_bins.R $1 $2 > $3

SCRIPTS_DIR=$1
BAM=$2
COUNT=$3
OUTPUT=$4

Rscript $SCRIPTS_DIR/make_bins.R $BAM $COUNT > $OUTPUT
