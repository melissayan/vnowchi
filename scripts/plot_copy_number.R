#!/usr/bin/Rscript

#--------------------------------------------------------------------------
# Plot Copy Numbers for 1 Sample
#--------------------------------------------------------------------------
# Plots the ratio of observed over expected read counts in each of the 
# supplied bins for the whole sample as well as individual chromosomes
# are generated and placed in a directory with the base name of the 
# file plus _plots
#
# Input:
#   - <sampleID>.cn.txt  Copy number for a given sample
# Output:
#   - <sampleID>_<chrName|all>.png  Plots for each chromosome and the 
#                                   whole genome of a sample
# 
# Usage: Rscript plot_copy_number.R <bin_count_file.cn>
#
# Example: Rscript plot_copy_number.R \
#          'D:/Box Sync/Rhesus_Embryos/Oocyte_CN/500/rh150409-1-b1-B1_S1.cn'
#
#
# Modified by Melissa Yan to accommodate different genomes and plot continuous CN black bars
# Original script by Nathan Lazar:
#   - Original script:
#       Project Title: Oocyte_CN
#       Script Title: plot_copy_number.R
#       Author: Nathan Lazar (nathan dot lazar at gmail dot com)
#       Date: March 7, 2017
#       Availability: https://github.com/nathanlazar/Oocyte_CN/blob/master/plot_copy_number.R
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

out.dir <- sub('.cn.txt', '', cn.file)

name <- sub('.cn.txt', '', cn.file)
split.name <- strsplit(name, split='/', fixed=T)[[1]]
name <- split.name[length(split.name)]

if(grepl('500', cn.file)) bin.size=500
if(grepl('1000', cn.file)) bin.size=1000
if(grepl('2000', cn.file)) bin.size=2000
if(grepl('4000', cn.file)) bin.size=4000

counts <- read.table(cn.file, header=T, stringsAsFactors=F)
counts <- counts[complete.cases(counts$chr),]

# Loop to make the counts.seg object. There's probably a way better way to do this
# with dplyr
counts.seg <- data.frame(chr=counts$chr[1], cn=counts$cn[1], 
                         start=counts$idx[1], end=counts$idx[1],
                         stringsAsFactors = F)
j <- 1
for(i in 2:nrow(counts)) {
  chr <- counts$chr[i]
  start <- counts$idx[i]
  end <- counts$idx[i]
  cn <- counts$cn[i]
  if(chr != counts.seg$chr[j]) {  # If new chrom
    j <- j+1
    counts.seg <- rbind(counts.seg, c(chr, cn, start, end))
  } else {
    if(cn == counts.seg$cn[j]) { # If not a new bin
      counts.seg$end[j] <- counts$idx[i]
    } else { # If a new bin
      j <- j+1
      counts.seg <- rbind(counts.seg, c(chr, cn, start, end))
    }
  }
}
counts$chr <- factor(counts$chr, levels=c(chr.Keep))
counts.seg$chr <- factor(counts.seg$chr, levels=levels(counts$chr))
counts.seg$cn <- as.numeric(counts.seg$cn)
counts.seg$start <- as.numeric(counts.seg$start)
counts.seg$end <- as.numeric(counts.seg$end)
for (i in 2:(nrow(counts.seg)-1)) {
  if (!is.na(counts.seg$start[i])){
    counts.seg$start[i] <- counts.seg$end[i-1]
    counts.seg$end[i] <- counts.seg$start[i+1]
  } 
}

# Prep for plotting
bins.per.chrom <- counts %>% tbl_df %>% group_by(chr) %>% summarise(len=length(chr)) 
bins.per.chrom$mid <- cumsum(bins.per.chrom$len) - (bins.per.chrom$len/2)
n.chr <- length(unique(counts$chr))
counts$col <- 2
counts$col[counts$chr %in% unique(counts$chr)[seq(1,n.chr,2)]] <- 1
counts$col <- as.factor(counts$col)
ylim <- c(0, min(max(counts$ratio), 6))

# Plot whole genome
png(file=paste0(out.dir, '_all.png'), width = 480*2)
p <- ggplot(counts, aes(x=idx, y=ratio)) +
  geom_point(aes(colour=col)) +
  theme(legend.position="none") +
  labs(y='Copy Number', x='Chromosome') +
  scale_x_continuous(breaks=cumsum(bins.per.chrom$len), 
                     minor_breaks=NULL, labels=rep("", nrow(bins.per.chrom))) +
  scale_y_continuous(limits=ylim, breaks=0:6,
                     minor_breaks=NULL) +
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
  labs(title=name) +
  geom_segment(data=counts.seg, aes(x=start, xend=end, y=cn, yend=cn), size=2)

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
dev.off()

# Plot each chromosome
for(chr in unique(counts$chr)) {
  chr.counts <- counts[counts$chr==chr,]
  chr.counts$idx <- chr.counts$bin
  chr.counts.seg <- counts.seg
  chr.counts.seg <- counts.seg[counts.seg$chr==chr,]
  chr.start <- min(chr.counts.seg$start)
  chr.counts.seg$start <- chr.counts.seg$start - chr.start
  chr.counts.seg$end <- chr.counts.seg$end - chr.start
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
  
  png(file=paste0(out.dir, '_', chr, '.png'), width = 480*2)
  p <- ggplot(chr.counts, aes(x=idx, y=ratio)) +
    geom_point(aes(colour=col)) +
    theme(legend.position="none") +
    labs(y='Copy Number', x="10Mb intervals") +
    scale_x_continuous(breaks=c(0,cumsum(mb10.bins$mb.bins)[-nrow(mb10.bins)]), 
                       labels=mb10.bins$mb) +
    scale_y_continuous(limits=ylim, breaks=0:6,
                       minor_breaks=NULL) +
    scale_colour_manual(values=cols) +
    theme(plot.title=element_text(size=24, face="bold", hjust=.5),
          axis.title=element_text(size=20),
          axis.text.x=element_text(size=20),
          axis.ticks.x=element_blank(), 
          axis.text.y=element_text(size=20)) +
    labs(title=paste(name, 'chromosome', chr.num)) +
    geom_segment(data=chr.counts.seg, aes(x=start, xend=end, y=cn, yend=cn), size=2)
  print(p)
  dev.off()
}


