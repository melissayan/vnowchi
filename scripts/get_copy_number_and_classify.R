#!/usr/bin/Rscript

#--------------------------------------------------------------------------
# Get Copy Numbers for Plots and Classify Sample Ploidy and Sex based on 
# Copy Number Variations
#--------------------------------------------------------------------------
# Uses CBS algorithm with package DNAcopy to determine regions of putative
# gain or putative loss and to estimate copy numbers given ratios
# calculated with count_bin_reads.R. Copy number estimates are 
# an added column in the given matrix of counts and ratios. The HMM
# algorithm with package HMMcopy is used to determine regions of putative
# gain or loss for comparison with CBS to identify copy number variant
# (CNVs) calls based on Knouse et al., 2016 and Knouse et al., 2017.
#
# Input: 
#   - <sampleID>.counts       counts file
#   - <genome>_chr_sizes.txt  chromosome size file
#   - VNOWCHI_summary.txt     statistical summary table
# Output: 
#   - <sampleID>.cn.txt   estimate of copy numbers based on ratios file
#   - <sampleID>.all.txt  Grouped positions based on gains and losses before filtering
#   - <sampleID>.CNV.txt  CNV regions and categorize sample ploidy and sex 
#
# Usage: Rscript get_copy_number.R <sampleID.counts> <genome_chr_sizes.txt> \
#        <VNOWCHI_summary.txt>
#
# Example: Rscript get_copy_number.R \
#          'D:/Box Sync/Rhesus_Embryos/Oocyte_CN/500/rh150409-1-b1-B1_S1.counts' \
#          ~~~~~~~~~ FINISH EXAMPLE? 
#
#
# Written by Melissa Yan to accommodate different genomes, include HMM, CBS 
# and HMM intersection, and classify sample ploidy and sex based on CNVs.
# Built on original script from get_copy_number.R by Nathan Lazar and using
# HMM example from SCI-seq_CBS_HMM_calling.pl by Kristof Torkenczy:
#   - Original script:
#       Project Title: Oocyte_CN
#       Script Title: get_copy_number.R
#       Author: Nathan Lazar (nathan dot lazar at gmail dot com)
#       Date: March 1, 2017
#       Availability: https://github.com/nathanlazar/Oocyte_CN/blob/master/get_copy_number.R 
#   - HMM Example:
#       Project Title: SCI-DNA_pipeline
#       Script Title: SCI-seq_CBS_HMM_calling.pl
#       Author: Kristof Torkenczy
#       Date: December 18, 2017
#       Code version: version 052616
#--------------------------------------------------------------------------

library(dplyr)    # General data frame manipulation
library(DNAcopy)  # Implements circular binary segmentation (CBS)
library(HMMcopy)  # Implements hidden markov model (HMM)
library(IRanges)  # Required for HMM
library(data.table)

set.seed(123) # Ensures CBS permutations are the same everytime

#--------------------------------------------------------------------------
# Change chromosomes names and size file to match the genome used
# Ex. RheMac8 has chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11
#                 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
#                 chrX chrY chrM
#     so copy following below this box:
#       # Chromosome names and sizes in RheMac8
#       chr.Keep <- c(paste0('chr',1:20), 'chrX', 'chrY', 'chrM')
#       chrX <- 'chrX'
#       chrY <- 'chrY'
#       chrM <- 'chrM'
#       # Chromosome sizes based on https://www.ncbi.nlm.nih.gov/genome/215
#       size.file <- "Mmul8.0.1_chr_sizes.txt"
#       
# Ex. UMD3.1 has 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
#                23 24 25 26 27 28 29 X Y MT
#     so copy following below this box:
#       # Chromosome names in UMD3.1
#       chrKeep <- c(paste0('',1:29), 'X', 'Y', 'MT')
#       chrX <- 'X'
#       chrY <- 'Y'
#       chrM <- 'MT'
#       # Chromosome sizes based on https://www.ncbi.nlm.nih.gov/genome/82
#       size.file <- "UMD3.1.1_chr_sizes.txt"
#     
#--------------------------------------------------------------------------
chr.Keep <- 
chrX <- 
chrM <- 
chrY <- 
size.file <- 


#--------------------------------------------------------------------------
# Functions:
#   - get_seg()
#   - recalibrate_CBS()
#--------------------------------------------------------------------------
# Get CBS segmentation for each chromosome
get_seg <- function(counts, base.cn) {
  # Use circular binary segmentation to call copy number states 
  counts.cna <- CNA(genomdat=counts$ratio * base.cn, 
                    chrom=counts$chr, 
                    maploc=1:nrow(counts), 
                    presorted=T)
  
  # This smoothes outliers. If a bin's ratio is more than 3 std. deviations
  # from it's nearest neighbor it will be shrunk to be 2 s.d. from the median.
  # The s.d. is taken from the distribution of the nearest 20 bins
  # Trim trims the most outlying 5% of data
  counts.smooth <- smooth.CNA(counts.cna, 
                              smooth.region=10, 
                              outlier.SD.scale=4,
                              smooth.SD.scale=2,
                              trim=0.025)
  
  # Run the circular binary segmentation.
  # Alpha is the significance level to accept change points.
  # nperm is the number of permutations run.
  # min.width is the smallest number of bins for a copy number segment
  # Prune undoes splits when the increase in sum of squares for unsplit 
  # is less than 5%
  counts.seg <- segment(counts.smooth, weights=rep(1,nrow(counts)), 
                        alpha = 1e-6, min.width=5,
                        nperm=1e6, undo.splits='sdundo',
                        undo.SD=.25, verbose=1)
 
#  counts.seg <- segment(counts.smooth, weights=rep(1,nrow(counts)), 
#                        alpha=0.0001, min.width=2,
#                        nperm=10000, undo.splits='sdundo',
#                        undo.SD=3, verbose=1)
 
  # Categorize copy numbers by putative "gain", putative "loss", or 
  # "exclude" based on segment mean (Knouse et al., 2016 and
  # Knouse et., al 2017):
  #  - gain: segment mean > 1.32+1
  #  - loss: segment mean < 0.60+1
  #  - exclude: segment mean <= 1.32+1 and segment mean >= 0.60+1
  # *note: need to add 1 because CBS is centered around 1; Adding 1 seems to 
  # result in a similar CBS gain/loss cutoffs used in McConnell et al 2013. 
  # The gains/loss can also be obtained by comparing CBS median copy number
  # by at least 2 median absolute deviation of predicted copy number values
  # in autosomal genomic windows (McConnell et al 2013) 
  counts.seg$output$round <- round(counts.seg$output$seg.mean)
  counts.seg$output$GL <- ifelse(counts.seg$output$seg.mean > 2.32, "gain", 
                                 ifelse(counts.seg$output$seg.mean < 1.60, "loss", "exclude"))
  counts.seg
}

# Recalibrate normalization based on the CBS segmentation counts 
recalibrate_CBS <- function(counts, counts.seg, bin.size=200, base.cn=2, verbose=T) {
  # Assign putative copy number to each window
  copynumber <- mapply(rep, counts.seg$round, counts.seg$num.mark)
  counts$cn <- as.numeric(unlist(copynumber, use.names=TRUE))
  
  if(verbose) {
    print('Putative copy number means before correction:')
    counts %>% tbl_df %>% group_by(cn) %>%
      summarise(mean.ratio=mean(ratio)) %>% print
  }
  
  # Correct for the differing amount of DNA in aneuploid samples
  control.tot <- sum(counts$count!=0) * bin.size * base.cn
  samp.tot <- sum(counts$cn[counts$count!=0])* bin.size/base.cn * 2
  # If all bins are called as copy number 0, don't correct
  if(samp.tot == 0) samp.tot <- control.tot
  counts$ratio <- counts$ratio * (samp.tot/control.tot) * base.cn
  
  if(verbose) {
    print('Putative copy number means after correction:')
    counts %>% tbl_df %>% group_by(cn) %>%
      summarise(mean.ratio=mean(ratio)) %>% print
  }
  print('recalibrate CBS done ~~~~~~~~')
  counts
}




#-------------------------------------------------------------------------
# Main script begins here
#-------------------------------------------------------------------------
# Load data
args <- commandArgs(trailingOnly = TRUE)
count.file <- args[1]
size.file <- args[2]
map.file <- args[3]

if(grepl('500', count.file)) bin.size=500
if(grepl('1000', count.file)) bin.size=1000
if(grepl('2000', count.file)) bin.size=2000
if(grepl('4000', count.file)) bin.size=4000
counts <- read.table(count.file, header=T, stringsAsFactors=F)
counts$chr <- factor(counts$chr, levels=c(chr.Keep), ordered=TRUE)

sizes <- read.table(size.file, header=T, sep="\t", stringsAsFactors=F)
sizes$chr <- factor(sizes$chr, levels=chr.Keep, ordered=TRUE)

map <- read.table(map.file, header=T, sep="\t", stringsAsFactors=F)

# Output file names (<sampleName>.cn.txt and <sampleName>.CNV.txt) and sample name
out <- sub('.counts', '.cn.txt', count.file)
out.cnv <- sub('.counts', '.CNV.txt', count.file)
name <- sub('.counts', '', count.file)
split.name <- strsplit(name, split='/', fixed=T)[[1]]
name <- split.name[length(split.name)]

cat("\n", name, "\n\n")

# Remove M and Y chroms
tot.reads <- sum(counts$count)
reads.rem <- sum(counts$count[counts$chr %in% c(chrM, chrY)])
per.rem <- reads.rem/tot.reads
exp.per.rem <- sum(counts$expected[counts$chr %in% c(chrM, chrY)])
print(sprintf('%.2f%% of reads removed that mapped to chromosomes %s and %s', 
              per.rem * 100, chrM, chrY))
adj.fact <- (1-exp.per.rem)/(1-per.rem)
print(sprintf('Expected %.2f%% of reads to map to chromosomes %s and %s so we correct by a factor of %.2f', 
              exp.per.rem * 100, chrM, chrY, adj.fact))
counts <- counts[!(counts$chr %in% c(chrM, chrY)),]
counts <- droplevels(counts)

# Adjust and ratios for removed reads
counts$per <- counts$per * adj.fact
counts$ratio <- counts$per/counts$expected

# Reorder chromosomes
counts <- with(counts, counts[order(chr),])

# Set bins w/ less than 2 reads to zero
read.cut <- 2
per.rem <- sum(counts$count[counts$count < read.cut])/sum(counts$count)
print(sprintf('%.2f%% of reads removed in bins with less than %s reads', 
              per.rem * 100, read.cut))
counts$count[counts$count < read.cut] <- 0
# Adjust expected counts removing bins w/ zero counts
counts$expected <- counts$expected * (1-per.rem/sum(counts$count))
counts$ratio[counts$count < read.cut] <- 0
counts$per[counts$count < read.cut] <- 0

# If most bins are empty we expect a base copy number of 1 
# otherwise base copy number is set to 2
if(mean(counts$count == 0) > .5) {
  base.cn <- 1
} else base.cn <- 2

# Store ratios for HMM calling
gc <- counts
gc <- gc[!is.na(gc$chr), ]


#-------------------------------------------------------------------------
# CBS calling
#-------------------------------------------------------------------------
# Get segmentation for each chromosome separately
chroms <- levels(counts$chr)
seg.list <- list()
for(chrom in chroms) {
  count.chr <- filter(counts, chr==chrom)
  seg.list[[chrom]] <- get_seg(count.chr, base.cn)$output
}
counts.seg <- seg.list[[1]]
for(i in 2:length(seg.list)){
  counts.seg <- rbind(counts.seg, seg.list[[i]])
}

# Recalibrate counts given the putative copy numbers
counts <- recalibrate_CBS(counts, counts.seg, bin.size, base.cn, verbose=T)

# Renumber the segmentation x-values
for(i in 2:nrow(counts.seg)) {
  counts.seg$loc.start[i] <- counts.seg$loc.end[i-1]+1
  counts.seg$loc.end[i] <- counts.seg$loc.start[i] + counts.seg$num.mark[i] - 1
}

# Add index for whole genome and each chrom
counts$idx <- 1:nrow(counts)
bins.per.chrom <- counts %>% tbl_df %>% group_by(chr) %>% summarise(len=length(chr)) 
counts$bin <- unlist(sapply(bins.per.chrom$len, function(x) 1:x))

# Add mean, gains/losses, and CBS copy number into counts
counts$cn <- 0
for(i in 1:nrow(counts.seg)){
  counts$CBS.mean[counts$idx >= counts.seg$loc.start[i] & 
                    counts$idx <= counts.seg$loc.end[i]] <- counts.seg$seg.mean[i]
  
  counts$CBS.GL[counts$idx >= counts.seg$loc.start[i] & 
                  counts$idx <= counts.seg$loc.end[i]] <- counts.seg$GL[i]
}
for(i in 1:nrow(counts)){
  # CBS cn for gain
  if (counts$CBS.mean[i] > 2.32) {
    if (counts$CBS.mean[i] < 2.5) {
      counts$cn[i] <- ceiling(counts$CBS.mean[i])
    } else {
      counts$cn[i] <- round(counts$CBS.mean[i])
    }
  # CBS cn for loss
  } else if (counts$CBS.mean[i] < 1.60) {
    if (counts$CBS.mean[i] > 1.50) {
      counts$cn[i] <- floor(counts$CBS.mean[i])
    } else {
      counts$cn[i] <- round(counts$CBS.mean[i])
    }
  } else {
    counts$cn[i] <- 2
  } 
}
counts$CBS.cn <- counts$cn


#-------------------------------------------------------------------------
# HMM calling
#-------------------------------------------------------------------------
# Create a dataframe with IRanges for HMM calling.  Make sure the start 
# position is smaller than the end position (no negative widths are 
# allowed for HMMsegment). Use psuedo fixed width coordinates of 
# width=1.
temp.gc <- transform(gc, start=pmin(gc$start, gc$end), end=pmax(gc$start, gc$end))
temp.gc$ratio[temp.gc$ratio==0] <- 0.00001
n <- length(gc$start)
chromwin <- c()
for(val in levels(gc$chr)) {
  chrlen <- length(which(temp.gc$chr==val))
  chromwin <- append(chromwin,c(1:chrlen))
}
gcalt <- IRanges(start=chromwin, width=rep.int(1,times=n))
gc.rangedwindow <- RangedData(ranges=gcalt, 
                              space=temp.gc$chr, 
                              per=gc$per, 
                              expected=gc$expected, 
                              ratio=gc$ratio, 
                              copy=log2(temp.gc$ratio))

# Optimize parameters. Change e parameter to 0.995 (Knouse et al.,2016)
parameters <- HMMsegment(gc.rangedwindow, getparam = TRUE)
parameters$e <- rep(0.995,length(parameters$e))
parameters

# Do HMM segmentation
segmented_copy <- HMMsegment(gc.rangedwindow, param=parameters)

# Categorize copy numbers by putative "gain", putative "loss", or 
# "exclude" based on median log2 ratios for HMM (Knouse et al.,
#  2016 and Knouse et., al 2017):
#  - gain: median log2 ratio > 0.4
#  - loss: median log2 ratio < -0.35
#  - exclude: median log2 ratio <= 0.4 and median log2 ratio >= -0.35
n <- length(segmented_copy$segs$chr)
gc.HMM <- data.frame(segmented_copy$segs$chr, 
                     segmented_copy$segs$start, 
                     segmented_copy$segs$end, 
                     (segmented_copy$segs$end - segmented_copy$segs$start+1), 
                     segmented_copy$segs$median, 
                     as.numeric(as.character(segmented_copy$segs$state)), 
                     as.numeric(as.character(segmented_copy$segs$state))-1)
names(gc.HMM) <- c("chr","start","end","width","HMM.median", "HMM.state", "HMM.cn")
gc.HMM$HMM.GL <- ifelse(gc.HMM$HMM.median > 0.4, "gain", ifelse(gc.HMM$HMM.median < -0.35, "loss", "exclude"))

# Renumber the segmentation x-values
for(i in 2:nrow(gc.HMM)) {
  gc.HMM$start[i] <- gc.HMM$end[i-1]+1
  gc.HMM$end[i] <- gc.HMM$start[i] + gc.HMM$width[i] - 1
}

# Add index for whole genome and each chrom
gc$idx <- 1:nrow(gc)
bins.per.chrom <- gc %>% tbl_df %>% group_by(chr) %>% summarise(len=length(chr)) 
gc$bin <- unlist(sapply(bins.per.chrom$len, function(x) 1:x))

# Add median, state, copy number, and gains/losses into gc
gc$HMM.cn <- 0
for(i in 1:nrow(gc.HMM)){
  gc$HMM.median[gc$idx >= gc.HMM$start[i] & 
                  gc$idx <= gc.HMM$end[i]] <- gc.HMM$HMM.median[i]
  
  gc$HMM.state[gc$idx >= gc.HMM$start[i] & 
                 gc$idx <= gc.HMM$end[i]] <- gc.HMM$HMM.state[i]
  
  gc$HMM.GL[gc$idx >= gc.HMM$start[i] & 
              gc$idx <= gc.HMM$end[i]] <- gc.HMM$HMM.GL[i] 
  
  gc$HMM.cn[gc$idx >= gc.HMM$start[i] & 
              gc$idx <= gc.HMM$end[i]] <- gc.HMM$HMM.cn[i]
}


#-------------------------------------------------------------------------
# Combine CBS and HMM 
#-------------------------------------------------------------------------
# Combine CBS and HMM results
combined <-0
combined <- merge(counts, gc, 
                  by=c("chr", "start", "end", "count", "per", "expected", "idx", "bin"))
colnames(combined) <- c("chr", "start", "end", "count", "per", "expected", "idx", "bin",
                        "ratio", "cn", "CBS.mean", "CBS.GL", "CBS.cn",
                        "HMM.ratio", "HMM.cn", "HMM.median", "HMM.state", "HMM.GL")
combined$chr <- factor(combined$chr, levels=c(chr.Keep))
combined <- arrange(combined, chr, start)

# If sample has UniqueMappingPositions < 50000, then the sample ploidy status is "no WGA"
if (nrow(map[map$SampleName==name, ]) == 0){
  combined$map.flag <- "missing.info"
} else if (!is.null(map[map$SampleName==name, ]$UniqueMappingPositions) &&
    (map[map$SampleName==name, ]$UniqueMappingPositions < 50000)){
  combined$map.flag <- "no.WGA"
} else {
  combined$map.flag <- NA
}


#-------------------------------------------------------------------------
# Add CBS and HMM putative gains and losses combined result (this will
# exclude euploid results where CNV=2)
#-------------------------------------------------------------------------
cnv <- 0
cnv <- combined[ ,c("chr", "start", "end", "count", "per", "expected", "ratio", "cn", "idx", "bin",
                    "CBS.mean", "CBS.cn", "CBS.GL",
                    "HMM.ratio", "HMM.median", "HMM.cn", "HMM.GL", "map.flag")]
cnv$intersect.GL <- interaction(cnv$CBS.GL, cnv$HMM.GL)
cnv$CBS.vs.HMM <- ifelse((cnv$intersect.GL=="loss.loss" | cnv$intersect.GL=="gain.gain"), "agree", "disagree")

# Make sure the start position is smaller than the end position
# Sometimes the windows at the end of chromosomes will end with a smaller position than the start position
chrom.list<-sort(union(levels(cnv$chr), levels(sizes$chr)))

cnv<-inner_join(mutate(cnv, chr=factor(chr, levels=chrom.list)),mutate(sizes, chr=factor(chr, levels=chrom.list)),by="chr")
cnv$end <- ifelse(cnv$end < cnv$start, cnv$size.in.bases, cnv$end)
cnv$chr <- factor(cnv$chr, levels=c(chr.Keep), ordered=TRUE)
cnv <- arrange(cnv, chr, start)


for(i in 1:nrow(cnv)){
  if (cnv$CBS.vs.HMM[i] == "disagree") {
    counts$cn[i] <- 2
    cnv$cn[i] <- 2
  }
}

# Save cnv data.frame for all CBS and HMM results per window and for CBS CNV plots
write.table(cnv, file=out, quote=FALSE, sep="\t", row.names=FALSE)


#-------------------------------------------------------------------------
# Condense table based on genome coordinates and similarities in CBS and HMM gain/loss
# https://stackoverflow.com/questions/36634182/r-compare-rows-of-data-frame-to-combine-based-on-conditions
#-------------------------------------------------------------------------
cnv <- as.data.table(cnv)

#If sample is euploid, create empty table with 1 row of NA
if (nrow(cnv) == 0){
  cnv <- setNames(data.frame(matrix(ncol=10, nrow=1)),
           c("chr", "start", "end", "width", "CBS.GL", "HMM.GL",
             "intersect.GL", "CBS.vs.HMM", "map.flag"))
} else {
  cnv <- cnv[, .(start = min(start), end = max(end), CBS.GL=CBS.GL[1], HMM.GL=HMM.GL[1],
                 intersect.GL=intersect.GL[1], CBS.vs.HMM=CBS.vs.HMM[1], map.flag=map.flag[1]),
             by = .(chr, rleid(intersect.GL), rleid(CBS.vs.HMM))]
  cnv <- as.data.frame(cnv)
  cnv$width <- (cnv$end - cnv$start)
  cnv <- cnv[c("chr", "start", "end", "width", "CBS.GL", "HMM.GL",
               "intersect.GL", "CBS.vs.HMM", "map.flag")]
  cnv <- data.frame(sample.name=name, cnv)
  cnv$chr <- factor(cnv$chr, chr.Keep, ordered=TRUE)
  cnv <- arrange(cnv, chr, start)
}

#DO NOT USE THIS. Dplyr doesn't work as intended on the server.  Could plyr R package be interferring?
#cnv <- cnv %>% dplyr::group_by(chr, CBS.GL, HMM.GL, intersect.GL, CBS.vs.HMM, map.flag) %>%
#  do( as(as.data.frame(.), "GRanges") %>% 
#        GenomicRanges::reduce() %>% 
#        as.data.frame() %>% 
#        dplyr::select(start, end)
#  ) %>%
#  ungroup() %>%
#  arrange(chr, start)

#DO NOT USE THIS. Genomic position gets grouped out of order
#cnv <- cnv %>% group_by(chr, CBS.GL, HMM.GL, CBS.vs.HMM) %>%
#  summarise(start=as.integer(paste(min(unique(start)),collapse=", ")), end=as.integer(max(end)),
#            count.sum=sum(count), cn.list=paste(sort(unique(cn)),collapse=", "),avg.CBS.ratio=mean(ratio), avg.HMM.ratio=mean(HMM.ratio),
#            CBS.cn.list=paste(sort(unique(CBS.cn)),collapse=", "), CBS.mean.list=paste(sort(unique(CBS.mean)),collapse=", "),
#            HMM.cn.list=paste(sort(unique(HMM.cn)),collapse=", "), HMM.median.list=paste(sort(unique(HMM.median)),collapse=", "),
#            map.flag=paste(sort(unique(map.flag)),collapse=", "),intersect.GL=paste(sort(unique(intersect.GL)),collapse=", "))
#cnv$width <- (cnv$end - cnv$start)
#cnv <- cnv[c("chr", "start", "end", "width", "count.sum", "cn.list", 
#             "avg.CBS.ratio", "avg.HMM.ratio",
#             "CBS.mean.list", "CBS.cn.list", "CBS.GL",
#             "HMM.median.list", "HMM.cn.list", "HMM.GL", 
#             "map.flag","intersect.GL","CBS.vs.HMM")]


#-------------------------------------------------------------------------
# Classify sample's ploidy status and sex based on regions of agreeing
# putative gaines or putative losses
#-------------------------------------------------------------------------
if (nrow(cnv) != 0){
  # Store disagreeing regions for later use
  cnv.no <- filter(cnv, cnv$CBS.vs.HMM == "disagree")

  # Get only agreeing regions
  cnv <- filter(cnv, cnv$CBS.vs.HMM == "agree")
}

# Classify ploidy of sample
if (all(is.na(cnv$sample.name)) && nrow(cnv) == 0 ){
  cnv[1, c("sample.name")] <- name
  cnv[1, c("sample.ploidy")] <- "euploid"
  cnv[1, c("flag")] <- NA
} else if (!is.na(cnv[1,c("map.flag")]) &&
           length(which(cnv[1,c("map.flag")] == "no.WGA")) ==1){
  cnv$sample.ploidy <- cnv$map.flag
  cnv$flag <- "check-No.WGA"
} else if (!is.null(length(which(is.element(chroms, cnv$chr)))) &&
           length(which(is.element(chroms, cnv$chr))) >=1 &&
           length(which(is.element(chroms, cnv$chr))) <=4 ) {
  cnv$sample.ploidy <- "Aneuploid"
  cnv$flag <- NA
} else if (!is.null(length(which(is.element(chroms, cnv$chr)))) &&
           length(which(is.element(chroms, cnv$chr)) > 5)) {
  cnv$sample.ploidy <- "chaotic.aneuploid"
  cnv$flag <- NA
}

# Classify sex of sample
if (all(subset(counts,chr==chrX)$cn == 2)) {
  cnv$sample.sex <- "female"
} else if (all(subset(counts,chr==chrX)$cn == 1)) {
  cnv$sample.sex <- "male"
} else {
  cnv$sample.sex <- "unknown"
}


#-------------------------------------------------------------------------
#Combine agreeing and diagreeing regions
#-------------------------------------------------------------------------
#Attach ploidy, sex, and flag results to the non-cnv data
if (nrow(cnv.no) != 0 ){
  cnv.no$sample.ploidy <- unique(cnv$sample.ploidy)
  cnv.no$flag <- unique(cnv$flag)
  cnv.no$sample.sex <- unique(cnv$sample.sex)
}

#Combine all data together
cnv <- bind_rows(cnv, cnv.no)


#-------------------------------------------------------------------------
# Determine if the impact on chromosomal regions are whole or segmental
# and reclassify the sample if it's aneuploid based on chromosomal regions
#-------------------------------------------------------------------------
# Determine chromosomal region's impact range
if (length(which(cnv[1,c("sample.ploidy")] != "euploid")) == 1){
  cnv$chrom.affected <- cnv$chr
  chrom.list<-sort(union(levels(cnv$chr), levels(sizes$chr)))
  cnv<-inner_join(mutate(cnv, chr=factor(chr, levels=chrom.list)),mutate(sizes, chr=factor(chr, levels=chrom.list)),by="chr")
  cnv$whole.or.seg <- ifelse(cnv$width > cnv$size.in.bases*0.9, "Whole", ifelse(cnv$width > 15000000, "Large Segment", "Small Segment"))
  cnv$gain.or.loss <- ifelse(cnv$intersect.GL == "gain.gain", "Gain", ifelse(cnv$intersect.GL == "loss.loss", "Loss", "Exclude"))
} else {
  cnv$chrom.affected <- NA
  cnv$size.in.bases <- NA
  cnv$whole.or.seg <- NA
  cnv$gain.or.loss <- NA
}

#If the sample is Aneuploid, then further categorize it as segmental or whole
if (length(which(cnv[1,c("sample.ploidy")] == "Aneuploid")) == 1){
  if (length(subset(subset(cnv, CBS.vs.HMM == "agree"), whole.or.seg == "Whole")$whole.or.seg) > 0){
    cnv$sample.ploidy <- "aneuploid.whole"
  } else {
    cnv$sample.ploidy <- "aneuploid.seg"
  }
}

#-------------------------------------------------------------------------
# Save data of final CBS and HMM results with ploidy and sex status
#-------------------------------------------------------------------------
#Reorganize data by chromosome and start position
cnv$chr <- factor(cnv$chr, chr.Keep, ordered=TRUE)
cnv <- arrange(cnv, chr, start)

# Reorganize columns
cnv <- cnv[c("sample.name", "chr", "start", "end", "width", 
             "CBS.GL", "HMM.GL", "CBS.vs.HMM",
             "chrom.affected", "size.in.bases", "whole.or.seg",
             "gain.or.loss", "sample.ploidy", "sample.sex", "flag")]

# Save cnv data.frame for condensed CBS and HMM results with ploidy and sex status
write.table(cnv, file=out.cnv, quote=FALSE, sep="\t", row.names=FALSE)

#Logs R session information for reproducibility purposes
sessionInfo()

cat("\n", name, "\t", cnv[1,c("sample.ploidy")], "\t", cnv[1,c("sample.sex")], "\t", cnv[1,c("whole.or.seg")], "\n")


