#!/usr/bin/Rscript

#--------------------------------------------------------------------------
# Make Bins Using Fibroblasts
#--------------------------------------------------------------------------
# Make a bed file of variable sized bins. Each bin should contain 
# approximately the same number of reads. 
#
# Input:
#   - xx-euploid-1_all_merge.bam SE file or xx-euploid-1_all.q30.sort.bam PE file
#   - <reads per bin>   number of reads per bin   
# Output:
#   - xx-euploid-1_all_<binSize>.bins or xx-euploid-1_all_merge_<binSize>.bins
#
# Usage: make_bins.R <reads>.bam <number of reads per bin>
#
# Example: make_bins.R MAPPED_MASKED/xx-euploid_all.q30.rmdup.bam 1000
#
#
# Modified by Melissa Yan to accommodate different genomes.
# Original script by Nathan Lazar:
#   - Original script:
#       Project Title: Oocyte_CN
#       Script Title: make_bins.R
#       Author: Nathan Lazar (nathan dot lazar at gmail dot com)
#       Date: February 4, 2017
#       Availability: https://github.com/nathanlazar/Oocyte_CN/blob/master/make_bins.R
#--------------------------------------------------------------------------

library(GenomicAlignments)

set.seed(123)

#--------------------------------------------------------------------------
# Change unmapped contigs to match the genome used
# Ex. RheMac8 has unmapped contigs which all contain 'Un' (ex.chrUn_NW_015090757v1)
#    unmapped.contigs <- c('Un')
# Ex. UMD3.1 has unmapped contigs which all contain 'GJ' (ex.GJ060270.1)
#    unmapped.contigs <- c('GJ')
# Ex. if geneome has multiple different unmapped contigs with different patterns
#    unmapped.contigs <- ('diffContig1', 'diffContig2')
#--------------------------------------------------------------------------
unmapped.contigs <- 




#-------------------------------------------------------------------------
# Main script begins here
#-------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
bam.file <- args[1]
n <- as.numeric(args[2])

reads <- readGAlignments(bam.file)
mcols(reads)$mid <- (end(reads)+start(reads)) %/% 2

# Subset to just the major chromosomes (without unmapped contigs)
big.chr <- levels(seqnames(reads))
big.chr <- setdiff(big.chr, unique(grep(paste(unmapped.contigs,collapse="|"), big.chr, value=TRUE)))
reads <- reads[seqnames(reads) %in% big.chr]

# Split up reads by chromosome
read.list <- split(reads, as.character(seqnames(reads)))

# Print the file header
cat(paste0('chr', '\t', 'start', '\t', 'end', '\t', 'count', '\n'))

# For each chromosome, split reads into chunks with n reads (middles) in each
for(i in 1:length(read.list)) {
  
  # Sort reads by midpoint
  read.list[[i]] <- read.list[[i]][order(mcols(read.list[[i]])$mid)]

  end <- 0
  remaining <- read.list[[i]]

  while(length(remaining) > n) {
    start <- end+1
    bin <- remaining[1:n]
    end <- max(mcols(bin)$mid)
    in.bin <- which(mcols(remaining)$mid >= start &
                      mcols(remaining)$mid <= end)
    bin <- remaining[in.bin]
    cat(paste0(names(read.list)[[i]], '\t', start, '\t', end, '\t', length(bin), '\n'))
    remaining <- remaining[-in.bin]
  }

  start <- end+1
  end <- seqlengths(read.list)[[i]]
  bin <- remaining
  cat(paste0(names(read.list)[[i]], '\t', start, '\t', end, '\t', length(bin), '\n'))
}


