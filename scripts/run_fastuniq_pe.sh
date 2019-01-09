#!/bin/bash

F1=$1
PATH=$PATH/:$2
export PATH

F2=${F1/R1/R2}
name1=${F1/.fq.gz/}; name1=$(basename $name1)
name2=${name1/R1/R2}
gunzip $F1
gunzip $F2

# Running FastUniq
echo ${F1/.gz/}  > COLLAPSED_PE/$name1.txt
echo ${F2/.gz/} >> COLLAPSED_PE/$name1.txt
fastuniq -i COLLAPSED_PE/$name1.txt -o COLLAPSED_PE/$name1.uniq.fq -p COLLAPSED_PE/$name2.uniq.fq

# Clean up 
rm COLLAPSED_PE/$name1.txt
gzip ${F1/.gz/}
gzip ${F2/.gz/}
gzip COLLAPSED_PE/$name1.uniq.fq
gzip COLLAPSED_PE/$name2.uniq.fq
