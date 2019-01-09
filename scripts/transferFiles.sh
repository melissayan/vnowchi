#!/usr/bin/bash

#Change these 2 items:
#  FILE_DIR is the location of the files to transfer
#  OUTPUT_DIR is the destination of where the files will transfer to
FILE_DIR='/home/groups/hoolock/u1/chavez_data/brittany_data/CopyNumberPipeline/MAPPED_PE'
OUTPUT_DIR='/home/exacloud/lustre1/yanm/transfer'

#Change the file names 
#find "$FILE_DIR" -name '170725_Chavez.160217_2_D3_ExF1*.bam' -exec cp  {} "$OUTPUT_DIR"  \;

find "$FILE_DIR" -name '170725_Chavez.160217_2_D3_ExF1*.bam' -exec cp  {} "$OUTPUT_DIR"  \;



