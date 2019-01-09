#!/bin/bash

f=$1
PATH=$PATH/:$2
export PATH

name=${f/.fq.gz/}; name=$(basename $name)
gunzip -c $f | \
fastx_collapser -Q33 > COLLAPSED_SE/$name.uniq.fq
gzip COLLAPSED_SE/$name.uniq.fq

