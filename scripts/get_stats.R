#!/usr/bin/Rscript

# nathan dot lazar at gmail dot com

# Script to form summary statistics into a table

library(reshape2)

args<-commandArgs(TRUE)

if (length(args)<1 || length(args)>2) {
   stop("please provide proper file path to FASTQC/ and FASTQC_TRIM/ directories")
} else if (length(args)==1) {
   # For PIPELINE_VNOWC
   DIR <- args[1]
   
   # Fastqc
   fastqcPath <- file.path(DIR, "FASTQC", "summary.txt")
   fastqc <- read.table(fastqcPath, stringsAsFactors=F, 
                        sep='\t', comment.char='')

   # Trimmed fastqc
   trimPath <- file.path(DIR, "FASTQC_TRIM", "summary.txt")
   fastqc.trim <- read.table(trimPath, stringsAsFactors=F,
                              sep='\t', comment.char='')

   # Mapped
   mapPath1 <- file.path(DIR, "MAPPED_SE", "read_stats.txt")
   map.stats <- read.table(mapPath1, stringsAsFactors=F,
                           sep='\t', comment.char='')
   mapPath2 <- file.path(DIR, "MAPPED_PE", "read_stats.txt")
   map.stats2 <- read.table(mapPath2, stringsAsFactors=F,
                           sep='\t', comment.char='')
   map.stats <- rbind(map.stats, map.stats2)

} else if (length(args)==2) {
   #For FIBROBLASTS PIPELINE_bins
#   DIR <- file.path(args[1], args[2])
   DIR <- file.path(args[1])

   # Fastqc
   fastqcPath <- file.path(DIR, "FASTQC", "summary.txt")
   fastqc <- read.table(fastqcPath, stringsAsFactors=F, 
                        sep='\t', comment.char='')

   # Trimmed fastqc
   trimPath1 <- file.path(DIR, "FASTQC_TRIM_SE", "summary.txt")
   fastqc.trim <- read.table(trimPath1, stringsAsFactors=F,
                              sep='\t', comment.char='')
   trimPath2 <- file.path(DIR, "FASTQC_TRIM_PE", "summary.txt")
   fastqc.trim2 <- read.table(trimPath2, stringsAsFactors=F,
                              sep='\t', comment.char='')
   fastqc.trim <- rbind(fastqc.trim, fastqc.trim2)

   #Mapped
   mapPath1 <- file.path(DIR, "MAPPED_SE", "read_stats.txt")
   map.stats <- read.table(mapPath1, stringsAsFactors=F,
                           sep='\t', comment.char='')
   mapPath2 <- file.path(DIR, "MAPPED_PE", "read_stats.txt")
   map.stats2 <- read.table(mapPath2, stringsAsFactors=F,
                           sep='\t', comment.char='')
   map.stats <- rbind(map.stats, map.stats2)
}

# Get stats of fastqc data ----------------------------------------------- 

# Clean up file names
names(fastqc)[1] <- 'file'
fastqc$file <- sub('.fq_fastqc/fastqc_data.txt:', '', fastqc$file)
fastqc$file <- sub('_R1_001', '', fastqc$file)
fastqc$file <- sub('_R2_001', '', fastqc$file)
fastqc$file <- sub('.R1.', '', fastqc$file)
fastqc$file <- sub('.R2.', '', fastqc$file)
fastqc <- unique(fastqc)

# Reshape results
fastqc$type <- ''
fastqc$type[grepl('otal Sequences', fastqc$file)] <- 'Total'
fastqc$file <- sub('Total Sequences', '', fastqc$file)
fastqc$file <- sub('otal Sequences', '', fastqc$file)
fastqc$type[grepl('equence length', fastqc$file)] <- 'length'
fastqc$file <- sub('Sequence length', '', fastqc$file)
fastqc$file <- sub('equence length', '', fastqc$file)
fastqc$type[grepl('Total Deduplicated Percentage', fastqc$file)] <- 'fastqc.dedup.per'
fastqc$file <- sub('#Total Deduplicated Percentage', '', fastqc$file)
fastqc$file <- sub('Total Deduplicated Percentage', '', fastqc$file)
fastqc$type[grepl('GC', fastqc$file)] <- 'per.GC'
fastqc$file <- sub('%GC', '', fastqc$file)
fastqc$file <- sub('GC', '', fastqc$file)

stats <- dcast(file~type, data=fastqc, value.var='V2', fun.aggregate=function(x) x[1])


# Get stats of trimmed fastqc data --------------------------------------- 

# Clean up file names
names(fastqc.trim)[1] <- 'file'
fastqc.trim$file <- sub('.fq_fastqc/fastqc_data.txt:', '', fastqc.trim$file)
fastqc.trim$file <- sub('.trim', '', fastqc.trim$file)
fastqc.trim$file <- sub('.R1', '', fastqc.trim$file)
fastqc.trim$file <- sub('.R2', '', fastqc.trim$file)
fastqc.trim <- unique(fastqc.trim)
fastqc.trim$type <- ''

# Reshape results
fastqc.trim$type[grepl('otal Sequences', fastqc.trim$file)] <- 'Total.trimmed'
fastqc.trim$file <- sub('Total Sequences', '', fastqc.trim$file)
fastqc.trim$file <- sub('otal Sequences', '', fastqc.trim$file)
fastqc.trim$type[grepl('equence length', fastqc.trim$file)] <- 'length.trimmed'
fastqc.trim$file <- sub('Sequence length', '', fastqc.trim$file)
fastqc.trim$file <- sub('equence length', '', fastqc.trim$file)
fastqc.trim$type[grepl('#Total Deduplicated Percentage', fastqc.trim$file)] <- 'fastqc.dedup.per.trimmed'
fastqc.trim$file <- sub('#Total Deduplicated Percentage', '', fastqc.trim$file)
fastqc.trim$file <- sub('Total Deduplicated Percentage', '', fastqc.trim$file)
fastqc.trim$type[grepl('%GC', fastqc.trim$file)] <- 'per.GC.trimmed'
fastqc.trim$file <- sub('%GC', '', fastqc.trim$file)
fastqc.trim$file <- sub('GC', '', fastqc.trim$file)

stats2 <- dcast(file~type, data=fastqc.trim, value.var='V2', fun.aggregate=function(x) x[1])

all <- merge(stats, stats2, all=T)
all <- all[,c('file', 'Total', 'length', 'per.GC', 'fastqc.dedup.per', 
              'Total.trimmed', 'length.trimmed', 'per.GC.trimmed', 'fastqc.dedup.per.trimmed')]


# Get stats of mapping data ---------------------------------------------- 

# Find the categories of data in the long form file
splits <- which(grepl('###', map.stats$V1))
# print(map.stats[splits+1,])

total <- data.frame(file=sapply(strsplit(map.stats$V1[1:(splits[1]-1)], split=' '), function(x) x[1]),
  unique=sapply(strsplit(map.stats$V1[1:(splits[1]-1)], split=' '), function(x) x[2]))
mapped <- data.frame(file=sapply(strsplit(map.stats$V1[(splits[1]+2):(splits[2]-1)], split=' '), function(x) x[1]),
  mapped=sapply(strsplit(map.stats$V1[(splits[1]+2):(splits[2]-1)], split=' '), function(x) x[2]))
q30 <- data.frame(file=sapply(strsplit(map.stats$V1[(splits[2]+2):(splits[3]-1)], split=' '), function(x) x[1]),
  q30=sapply(strsplit(map.stats$V1[(splits[2]+2):(splits[3]-1)], split=' '), function(x) x[2]))
pos <- data.frame(file=sapply(strsplit(map.stats$V1[(splits[3]+2):nrow(map.stats)], split='  '), function(x) x[1]),
  pos=sapply(strsplit(map.stats$V1[(splits[3]+2):nrow(map.stats)], split='  '), function(x) x[2]))
m1 <- merge(total, mapped, all=T)
m2 <- merge(m1, q30, all=T)
m3 <- merge(m2, pos, all=T)

final <- merge(all, m3, all=T)

output <- file.path(DIR, "mapping_stats.txt")
write.table(final, file=output, quote=F, sep='\t', row.names=F)
