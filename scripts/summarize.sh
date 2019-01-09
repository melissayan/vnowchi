#!/usr/bin/bash

#--------------------------------------------------------------------------
# Save Mapped Read Statistics in Separate Files
#--------------------------------------------------------------------------
# Save the following in individual files:
#   - Total mapped reads
#   - Passing Q30 filtered reads
#   - Uniquely mapped reads
# These files will be used in generating summary read stats for:
#   - VNOWCHI_Summary.txt
#--------------------------------------------------------------------------
WORKING_DIR=$1

# Get mapped summary read stats for single-ended
cd $WORKING_DIR/MAPPED_SE

grep  Uniquely *.mapped_reads | sed 's/.mapped_reads:Uniquely mapped reads://; s/\.fastq.gz//; s/\.R1//' | sort | uniq | awk '{print $1 "\t" $2}' > 1_mapped.txt

grep Passing *.mapped_reads | sed 's/.mapped_reads:Passing Q 30 filter://; s/\.fastq.gz//; s/\.R1//' | sort | uniq | awk '{print $1 "\t" $2}' > 1_passQ30.txt

grep -A1 Total *.mapped_reads | awk '((NR-2)%3==0) {print $1,$2}' | sed -r 's/.mapped_reads-[0-9]+/ /; s/\.fastq.gz//; s/\.R1//' | sort | uniq | awk '{print $1 "\t" $2}' > 1_uniqueMappingPositions.txt

awk 'BEGIN{print "SampleName\tMapped\tPassQ30\tUniqueMappingPositions"} {a[$1]=a[$1]$2"\t"} END{for (i in a){print i "\t" a[i]}}' 1_mapped.txt 1_passQ30.txt 1_uniqueMappingPositions.txt > 1_temp.txt

(head -n 1 1_temp.txt && tail -n +2 1_temp.txt | sort) > summary_SE.txt


# Get mapped summary read stats for paired-end samples
cd $WORKING_DIR/MAPPED_PE

grep  Uniquely *.mapped_reads | sed 's/.mapped_reads:Uniquely mapped reads://; s/\.fastq.gz//; s/\.R1//; s/\.trim//g; s/\.fq//g; s/\.gz//g' | sort | uniq | awk '{print $1 "\t" $2}' > 1_mapped.txt

grep Passing *.mapped_reads | sed 's/.mapped_reads:Passing Q 30 filter://; s/\.fastq.gz//; s/\.R1//; s/\.trim//g; s/\.fq//g; s/\.gz//g' | sort | uniq | awk '{print $1 "\t" $2}' > 1_passQ30.txt

grep -A1 Total *.mapped_reads | awk '((NR-2)%3==0) {print $1,$2}' | sed -r 's/.mapped_reads-[0-9]+/ /; s/\.fastq.gz//; s/\.R1//g; ; s/\.trim//g; s/\.fq//g; s/\.gz//g' | sort | uniq | awk '{print $1 "\t" $2}' > 1_uniqueMappingPositions.txt

awk 'BEGIN{print "SampleName\tMapped\tPassQ30\tUniqueMappingPositions"} {a[$1]=a[$1]$2"\t"} END{for (i in a){print i "\t" a[i]}}' 1_mapped.txt 1_passQ30.txt 1_uniqueMappingPositions.txt > 1_temp.txt

(head -n 1 1_temp.txt && tail -n +2 1_temp.txt | sort) > summary_PE.txt

cd ..
cat MAPPED_SE/summary_SE.txt > summary_MAPPED.txt
tail -n +2 MAPPED_PE/summary_PE.txt >> summary_MAPPED.txt

cd $WORKING_DIR


