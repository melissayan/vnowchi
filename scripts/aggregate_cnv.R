#!/usr/bin/Rscript

library(dplyr)    # General data frame manipulation
library(tibble)

set.seed(123) 

#--------------------------------------------------------------------------
# Summarize Sample and Embryo Results
#--------------------------------------------------------------------------
# Uses results from get_copy_number_and_classify.R and information about
# embryo and samples to create 2 summary tables:
#    1. all sample copy number variation calls by embryo
#    2. embryo ploidy and sex status 
#
# Input:
#   - <sampleID>.CNV.txt  CNV regions and categorize sample ploidy and sex 
#   - group_info_SE.txt   Single-end sample name, embryo, #blastomere info
#   - group_info_PE.txt   Paired-end sample name, embryo, #blastomere info
# Output:
#   - CNV_<SE|PE>_<bin>.sampleSummary.txt   Summary of all individual
#                                           samples with embryo and
#                                           blastomere data
#   - CNV_<SE|PE>_<bin>.embryoSummary.txt   Summary of embryos with ploidy
#                                           and sex status
#
# Written by Melissa Yan 
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
# Change chromosomes to match the genome used
# Ex. RheMac8 has chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11
#                 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
#                 chrX chrY chrM
#   chr.Keep <- c(paste0('chr',1:20), 'chrX', 'chrY', 'chrM')
# Ex. UMD3.1 has 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
#                23 24 25 26 27 28 29 X Y MT
#    chrKeep <- c(paste0('',1:29), 'X', 'Y', 'MT')
#--------------------------------------------------------------------------
chr.Keep <- 


#--------------------------------------------------------------------------
# Functions:
#   - add_empty_column()
#--------------------------------------------------------------------------
# Adds empty columns to dataframe if columns don't exist in dataframe
# From: https://stackoverflow.com/questions/45857787/adding-column-if-it-does-not-exist
add_empty_column <- function(df, colnames){
  not.included.cols <- colnames[!colnames %in% names(df)]
  if(length(not.included.cols) != 0 ) {
    df[not.included.cols] <- 0
  }
  df
}




#--------------------------------------------------------------------------
# Main script begins here
#-------------------------------------------------------------------------
# Load data
args <- commandArgs(trailingOnly = TRUE)
cnv.file <- args[1]
group.file <- args[2]

cnv <- read.table(cnv.file, header=T, sep="\t", stringsAsFactors=F)
group <- read.table(group.file, header=T, sep="\t", stringsAsFactors=F)

# Output file names
out.sample.summary <- sub('.txt', '.sampleSummary.txt', cnv.file)
out.embryo.summary <- sub('.txt', '.embryoSummary.txt', cnv.file)


#--------------------------------------------------------------------------
# Summary of individual samples' CNV calls with embryo information
#--------------------------------------------------------------------------
# Add embryo group information to sample cnv and reorganize columns
cnv <- cnv %>% inner_join(group, by=c("sample.name" = "sample.name"))
cnv <- cnv[ ,c("embryo.name", "blastomeres", 
               "sample.name", "se.pe", "chr", "start", "end", "width", 
               "CBS.GL", "HMM.GL", "CBS.vs.HMM", 
               "chrom.affected", "size.in.bases", "whole.or.seg", 
               "gain.or.loss", "sample.ploidy", "sample.sex", "flag")]

# Filter information to show only agreeing regions (cnvs) in aneuploid
# samples, euploid samples, and no.WGA samples 
cnv.filtered <- cnv[cnv$CBS.vs.HMM == "agree" | (cnv$sample.ploidy == "euploid" & is.na(cnv$chr)) | (cnv$sample.ploidy == "no.WGA"), ]

write.table(cnv.filtered, file=out.sample.summary, quote=FALSE, sep="\t", row.names=FALSE)


#--------------------------------------------------------------------------
# Summary of embryos' ploidy and sex statuses
#--------------------------------------------------------------------------
# Summarize cnv to 1 sample per row 
embryos <- cnv[!duplicated(cnv$sample.name), ]

# Get embryo group information
group.info <- summarize(group_by(group,embryo.name,blastomeres))

# Add embryo ploidy information
ploidy <- as.data.frame.matrix(with(embryos, table(embryo.name,sample.ploidy)))
ploidy.names <- c("aneuploid.seg","aneuploid.whole", "chaotic.aneuploid","euploid", "no.WGA")
ploidy <- add_empty_column(ploidy, ploidy.names)
ploidy <- tibble::rownames_to_column(ploidy, var="embryo.name")
ploidy <- merge(group.info, ploidy, by="embryo.name")
ploidy <- ploidy[ ,c("embryo.name", "blastomeres", "aneuploid.seg", "aneuploid.whole", "chaotic.aneuploid","euploid", "no.WGA")]
for (i in 1:nrow(ploidy)) {
  row <- ploidy[i, ] 
  if (row$no.WGA > 0 & ((row$blastomeres - row$no.WGA) == 0)) {
    ploidy[i, c("embryo.ploidy")] <- "no.WGA"
  } else if ((row$blastomeres) == row$euploid) {
    ploidy[i, c("embryo.ploidy")] <- "euploid"
  } else if (row$aneuploid.seg>=0 && row$aneuploid.whole>=0 && row$chaotic.aneuploid>=0 && row$euploid==0) {
    ploidy[i, c("embryo.ploidy")] <- "aneuploid"
  } else if (row$aneuploid.seg>=0 && row$aneuploid.whole>=0 && row$chaotic.aneuploid>=0 && row$euploid>0){
    ploidy[i, c("embryo.ploidy")] <- "mosaic"
  } else {
    ploidy[i, c("embryo.ploidy")] <- "unknown"
  }
}

# Add embryo sex information
sex <- as.data.frame.matrix(with(embryos, table(embryo.name,sample.sex)))
sex.names <- c("female", "male", "unknown")
sex <- add_empty_column(sex, sex.names)
sex <- tibble::rownames_to_column(sex, var="embryo.name")
sex <- merge(group.info, sex, by="embryo.name")
sex <- as.data.frame.matrix(with(embryos, table(embryo.name,sample.sex)))
sex.names <- c("female", "male", "unknown")
sex <- add_empty_column(sex, sex.names)
sex <- rownames_to_column(sex, var="embryo.name")
sex <- merge(group.info, sex, by="embryo.name")
for (i in 1:nrow(sex)) {
  row <- sex[i, ]
  if (row$female == row$blastomeres && row$female > 0) {
    sex[i, c("embryo.sex")] <- "female"
    sex[i, c("flag")] <- ""
  } else if (row$male == row$blastomeres && row$male > 0){
    sex[i, c("embryo.sex")] <- "male"
    sex[i, c("flag")] <- ""
  } else {
    sex[i, c("embryo.sex")] <- "unknown"
    sex[i, c("flag")] <- "check"
  }
}

# Embryos with ploidy and sex information
embryos <- merge(ploidy, sex, by="embryo.name")
embryos$embryo.name <- factor(embryos$embryo.name, unique(group$embryo.name), ordered=TRUE)
embryos <- dplyr::rename(embryos, blastomeres = blastomeres.x )
embryos <- embryos %>% select(embryo.name, blastomeres, 
                   aneuploid.seg, aneuploid.whole, chaotic.aneuploid, euploid, no.WGA, embryo.ploidy,
                   female, male, unknown, embryo.sex, flag)

write.table(embryos, file=out.embryo.summary, quote=FALSE, sep="\t", row.names=FALSE)


#Logs R session information for reproducibility purposes
sessionInfo()

