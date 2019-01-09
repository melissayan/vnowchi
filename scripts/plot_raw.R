#!/usr/bin/Rscript

#--------------------------------------------------------------------------
# Plot Read Counts
#--------------------------------------------------------------------------
# Plots the read counts in each of the supplied bins for the
# whole sample as well as individual chromosomes. png files
# are generated and placed in a directory with the base name of the 
# file plus _plots
#
# Input:
#   - <sampleID>.cn.txt  Copy number for a given sample
# Output:
#   - <sampleID>_<chrName|all>.png  Plots for each chromosome and the 
#                                   whole genome of a sample
#
# Usage: Rscript plot_raw.R <bin_count_file.cn> <out_dir>
#
# Example: Rscript plot_raw.R \
#          'D:/Box Sync/Rhesus_Embryos/Oocyte_CN/4000/rh150409-1-b1-B1_S1.cn'
#
#
# Modified by Melissa Yan to accommodate different genomes.
# Original script by Nathan Lazar:
#   - Original script:
#       Project Title: Oocyte_CN
#       Script Title: plot_raw.R
#       Author: Nathan Lazar (nathan dot lazar at gmail dot com)
#       Date: March 7, 2017
#       Availability: https://github.com/nathanlazar/Oocyte_CN/blob/master/plot_raw.R
#--------------------------------------------------------------------------

library(dplyr)    # data frame manipulation
library(ggplot2)  # Used for plotting
library(grid)     # Used for plotting 

#-------------------------------------------------------------------------
# Change chromosomes to match the genome used
# Ex. RheMac8 has chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11
#                 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
#                 chrX chrY chrM
#    chr.Keep <- c(paste0('chr',1:20), 'chrX')
# Ex. UMD3.1 has 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
#                23 24 25 26 27 28 29 X Y MT
#    chrKeep <- c(paste0('',1:29), 'X')
#-------------------------------------------------------------------------
chr.Keep <- 




#-------------------------------------------------------------------------
# Main script begins here
#-------------------------------------------------------------------------
# Load data
args <- commandArgs(trailingOnly = TRUE)
cn.file <- args[1]
out.dir <- sub('.cn.txt', '_raw_plots', cn.file)

name <- sub('.cn.txt', '', cn.file)
split.name <- strsplit(name, split='/', fixed=T)[[1]]
name <- split.name[length(split.name)]

counts <- read.table(cn.file, header=T, stringsAsFactors=F, sep="\t")
counts <- counts[complete.cases(counts$chr),]

# Make output directory if necessary 
if (!file.exists(out.dir))
  dir.create(file.path(out.dir))

# Create factors for chromosomes in genome
counts$chr <- factor(counts$chr, levels=c(chr.Keep))

# Prep for plotting
bins.per.chrom <- counts %>% tbl_df %>% group_by(chr) %>% summarise(len=length(chr)) 
bins.per.chrom$mid <- cumsum(bins.per.chrom$len) - (bins.per.chrom$len/2)
n.chr <- length(unique(counts$chr))
counts$col <- 2
counts$col[counts$chr %in% unique(counts$chr)[seq(1,n.chr,2)]] <- 1
counts$col <- as.factor(counts$col)
ylim <- c(0, max(counts$count))

# Plot whole genome
png(file=paste0(out.dir, '/', name, '_all_raw.png'), width = 480*2)
p <- ggplot(counts, aes(x=idx, y=count)) +
  geom_point(aes(colour=col)) +
  theme(legend.position="none") +
  labs(y='Reads in each bin', x='Chromosome') +
  scale_x_continuous(breaks=cumsum(bins.per.chrom$len), 
                     minor_breaks=NULL, labels=rep("", nrow(bins.per.chrom))) +
  scale_y_continuous(limits=ylim, minor_breaks=NULL) +
  scale_colour_manual(values=c('#4861B3', '#70B333')) +
  theme(plot.title=element_text(size=24, face="bold", hjust=.5),
        axis.title=element_text(size=20),
        axis.text.x=element_text(size=14),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_text(size=20)) +
  geom_text(data=bins.per.chrom, size=6, colour='#7F7F7F',
            aes(x=bins.per.chrom$mid, vjust=2.5,
                y=rep(0,nrow(bins.per.chrom)), 
                label=sub('chr', '', bins.per.chrom$chr))) +
  labs(title=name) 

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
dev.off()

for(chr in unique(counts$chr)) {
  chr.counts <- counts[counts$chr==chr,]
  chr.counts$idx <- chr.counts$bin
  chr.bins.per.chrom <- bins.per.chrom[bins.per.chrom$chr==chr,]
  chr.bins.per.chrom$mid <- chr.bins.per.chrom$len/2
  chr.num <- gsub("chr", "", chr)

  # Hack to make colors alternate
  if(is.na(chr)) {
    cols <- c('#70B333', '#4861B3')
  } else if(unique(chr.counts$col)==1) {
    cols <- c('#4861B3', '#70B333')
  } else cols <- c('#70B333', '#4861B3')
  
  # Make column showing breaks every 10 Mb
  chr.counts$mb <- chr.counts$end %/% 10000000
  mb10.bins <- chr.counts %>% tbl_df %>% group_by(mb) %>% summarise(mb.bins=length(idx))
  
  png(file=paste0(out.dir, '/', name, '_', chr, '_raw.png'), width = 480*2)
  p <- ggplot(chr.counts, aes(x=idx, y=count)) +
    geom_point(aes(colour=col)) +
    theme(legend.position="none") +
    labs(y='Reads in each bin', x="10Mb intervals") +
    scale_x_continuous(breaks=c(0,cumsum(mb10.bins$mb.bins)[-nrow(mb10.bins)]), 
                       labels=mb10.bins$mb) +
    scale_y_continuous(limits=ylim, minor_breaks=NULL) +
    scale_colour_manual(values=cols) +
    theme(plot.title=element_text(size=24, face="bold", hjust=.5),
          axis.title=element_text(size=20),
          axis.text.x=element_text(size=20),
          axis.ticks.x=element_blank(), 
          axis.text.y=element_text(size=20)) +
    labs(title=paste(name, 'chromosome', chr.num))
  print(p)
  dev.off()
}
