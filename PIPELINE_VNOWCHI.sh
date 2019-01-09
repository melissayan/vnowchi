#!/usr/bin/bash
#
#SBATCH --get-user-env
#SBATCH --partition=very_long_jobs
#SBATCH --cpus-per-task=24
#SBATCH --mem=48000
#SBATCH --time=43200
#SBATCH --verbose
#SBATCH --output=/results/SLURM_debug/VNOWCHI_%j.out
#SBATCH --error=/results/SLURM_debug/VNOWCHI_%j.out
#SBATCH --export=ALL  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

set -e
set -x


#----------------------------------------------------------------------------------------
# Pipeline for Variable Non-Overlapping Window CBS & HMM Intersect (VNOWCHI)
#----------------------------------------------------------------------------------------
# Steps to map reads and create copy-number estimates for each sample. We use known 
# euploid samples to make variable-sized bins with a fixed expected number of reads
# in each bin. Then for each sample we count the number of reads in each bin and
# compare to what would be expected given the total number of reads.
#
# Circular Binary Segmentation (CBS) using the R package DNAcopy determines putative 
# gain or putative loss and calls putative changes in copy number which are adjusted 
# to whole numbers. Hidden Markov Model (HMM) using the R package HMMcopy also 
# determines putative gains or losses. Copy number variations (CNVS) are identified
# when CBS and HMM regions both identify a putative gain or putative loss (Knouse et 
# al., 2016 and Knouse et al., 2017).
#
# A base copy number of 2 is assumed, but set to 1 if most bins are empty. 
#
#
# Software versions:
#   - FastQC: v0.10.1
#   - Trimmomatic v0.35  *note: .jar file needs to be in specific directory see trim_dir
#   - FASTX-Toolkit: v0.0.13
#   - BWA-MEM: 0.7.9a-r786
#   - SAMtools: 0.1.19-44428cd
#   - BEDtools: v2.25.0
#   - FastUniq: v1.1
#
# Final output files of interest summarizing CNV per embryo based on individual samples
# for SE or PE based on bins:
#   - CNV_{SE|PE}_${bins}.sampleSummary.txt  Summary of all individual samples with
#                                            condensed CBS/HMM results, sample ploidy
#                                            status, and sample sex status.
#   - CNV_{SE|PE}_${bins}.embryoSummary.txt  Summary of embryo's ploidy and sex status
#   - VNOWCHI_Summary.txt
#
# References:
#   - Knouse, K. A., Wu, J., & Amon, A. (2016). Assessment of megabase-scale
#     somatic copy number variation using single-cell sequencing. Genome 
#     Research, 26(3), 376–384. http://doi.org/10.1101/gr.198937.115
#   - Knouse, K. A., Wu, J., & Hendricks, A. (2017). Detection of Copy Number
#     Alterations Using Single Cell Sequencing. Journal of Visualized 
#     Experiments: JoVE, (120), 55143. Advance online publication.
#     http://doi.org/10.3791/55143
#   - McConnell, M. J., Lindberg, M. R., Brennand, K. J., Piper, J. C., Voet, T.,
#     Cowing-Zitron, C., Shumilina, S., Lasken, R. S., Vermeesch, J. R., Hall, I. M., &
#     Gage, F. H. (2013). Mosaic copy number variation in human neurons. Science (New 
#     York, N.Y.), 342(6158), 632-7. http://doi.org/10.1126/science.1243472
#
# Modified by Melissa Yan for SLURM, to accommodate different genomes, include HMM, CBS 
# and HMM intersection, summary stats, and classify embryo ploidy and sex status.
#
# Built on original script by Nathan Lazar and HMM example from Kristof Torkenczy.
#   - Original script:
#       Project Title: Oocyte_CN
#       Script Title: PIPELINE_VNOWC
#       Author: Nathan Lazar (nathan dot lazar at gmail dot com)
#       Date: September 19, 2017
#       Availability: https://github.com/nathanlazar/Oocyte_CN/blob/master/PIPELINE_VNOWC 
#   - HMM Example:
#       Project Title: SCI-DNA_pipeline
#       Script Title: SCI-seq_CBS_HMM_calling.pl
#       Author: Kristof Torkenczy
#       Date: December 18, 2017
#       Code version: version 052616
#----------------------------------------------------------------------------------------


WORKING_DIR=$(pwd)/results
SCRIPTS_DIR=$(pwd)/scripts

#change this to your username
USERNAME=

# Add paths to tools
PATH=$PATH/:/home/users/$USERNAME/bin
PATH=$PATH/:/home/users/$USERNAME/bin/bedtools2/bin
PATH=$PATH/:/home/users/$USERNAME/bin/FastQC
PATH=$PATH/:/home/users/$USERNAME/bin/fastx/bin
PATH=$PATH/:/home/users/$USERNAME/bin/FastUniq/source
export PATH
trim_dir=/home/users/$USERNAME/bin/Trimmomatic-0.35

# Set the values for the Y chromosome and the mitochondrial DNA
chrY=
chrMT=

# Set the reference genome location
REF=

# Set the file name containing chromosome name and size
CHR_SIZE_TABLE=

# Set the file names for single-end and paired-end samples containing: 
#   - sample name
#   - embryo name
#   - number of blastomeres
GROUP_INFO_SE_TABLE="$(pwd)/group_info_SE.txt"
GROUP_INFO_PE_TABLE="$(pwd)/group_info_PE.txt"

# Set the number of parallel threads to use
cores=1

# Set srun to allocate single core to each task
SRUN="srun --exclusive -N1 -n1"

#----------------------------------------------------------------------------------------
# 0. Construct bins using fibroblast samples
#   See: PIPELINE_bins
#
# 1. Collect all raw read files to process
#    *note: The directory 'fastq' should have all .fq.gz files in either SE or PE for 
#           single-end files or paired-end files
#
#----------------------------------------------------------------------------------------
cd $WORKING_DIR

<<"COMMENT"
#COMMENT
#----------------------------------------------------------------------------------------
# 2. Run FastQC on the raw reads
#----------------------------------------------------------------------------------------
if [ ! -d FASTQC ]; then mkdir FASTQC; fi
if [ ! -d FASTQC/SE ]; then mkdir FASTQC/SE; fi
if [ ! -d FASTQC/PE ]; then mkdir FASTQC/PE; fi

# Run FastQC on SE and PE FASTQ files
for f in `find $WORKING_DIR/fastq/SE -name "*.gz"`; do
  $SRUN fastqc --outdir FASTQC/SE --threads $cores $f &
  #fastqc --outdir FASTQC/SE --threads $cores $f &
done
for f in `find $WORKING_DIR/fastq/PE -name "*.gz"`; do
  $SRUN fastqc --outdir FASTQC/PE --threads $cores $f &
  #fastqc --outdir FASTQC/PE --threads $cores $f &
done
wait

# Write FastQC results for each sample before trimming to summary_FASTQC.txt:
#    - Raw reads
#    - Read length 
#    - %GC
cd FASTQC/SE
for f in `ls *.zip`; do unzip $f; done
grep "Total Sequences" */fastqc_data.txt > 1_reads.txt
grep "Sequence length" */fastqc_data.txt > 1_readLength.txt
grep "^%GC" */fastqc_data.txt > 1_GC.txt
sed -i 's/ /\_/g; s/\_fastqc\/fastqc_data.txt\:/\t/g; s/\_R1//g; s/\_001//g; s/_L001//g; s/\.R1//g' *.txt
awk 'BEGIN{print "SampleName\tSeOrPe\tRawReads\tReadLength\t%GC"} {a[$1]=a[$1]$3"\t"} END{for (i in a){print i "\tSE\t" a[i]}}' 1_reads.txt 1_readLength.txt 1_GC.txt > 1_temp.txt
(head -n 1 1_temp.txt && tail -n +2 1_temp.txt | sort) > summary.txt
cd ../PE
for f in `ls *.zip`; do unzip $f; done
grep "Total Sequences" */fastqc_data.txt > 1_reads.txt
grep "Sequence length" */fastqc_data.txt > 1_readLength.txt
grep "^%GC" */fastqc_data.txt > 1_GC.txt
sed -i 's/ /\_/g; s/\_fastqc\/fastqc_data.txt\:/\t/g; s/\_R1//g; s/\_001//g; s/\.R1//g' *.txt
awk '{sum+=$3} NR%2==1{printf $1 "\t"} NR%2==0{print sum; sum=0;}' 1_reads.txt > 1_reads2.txt
awk 'NR%2==1{print $1 "\t" $3}' 1_readLength.txt > 1_readLength2.txt
awk 'NR%2==1{print $1 "\t" $3}' 1_GC.txt > 1_GC2.txt
awk 'BEGIN{print "SampleName\tSeOrPe\tRawReads\tReadLength\t%GC"} {a[$1]=a[$1]$2"\t"} END{for (i in a){print i "\tPE\t" a[i]}}' 1_reads2.txt 1_readLength2.txt 1_GC2.txt > 1_temp.txt
(head -n 1 1_temp.txt && tail -n +2 1_temp.txt | sort) > summary.txt
cd ..
cat SE/summary.txt > $WORKING_DIR/summary_FASTQC.txt
tail -n +2 PE/summary.txt >> $WORKING_DIR/summary_FASTQC.txt
rm $WORKING_DIR/FASTQC/*/1_*.txt
rm $WORKING_DIR/FASTQC/*/*.html
rm -r $WORKING_DIR/FASTQC/*/*fastqc
cd ..

echo ------------------- done FASTQC on raw reads;echo `date +"%Y-%m-%d %T"`;echo;echo


#----------------------------------------------------------------------------------------
# 3. Trim low quality bases w/ Trimmomatic
#----------------------------------------------------------------------------------------
if [ ! -d TRIM_SE ]; then mkdir TRIM_SE; fi
if [ ! -d TRIM_PE ]; then mkdir TRIM_PE; fi

# Trim single end reads
for F1 in `ls fastq/SE/*`; do
#  name=${F1/_R1_001.fastq.gz/}; name=${name/.R1.fastq.gz}; name=$(basename $name)
#  name=${F1/_R1_001.fq.gz/}; name=${name/.R1.fq.gz}; name=$(basename $name)
  name=${F1/.R1.gz/}; name=${name/.R1.fq.gz}; name=$(basename $name)
  $SRUN java -jar $trim_dir/trimmomatic-0.35.jar \
    SE -threads $cores -phred33 $F1 \
    TRIM_SE/$name.trim.fq.gz  \
    ILLUMINACLIP:$trim_dir/adapters/TruSeq3-SE.fa:2:30:10 LEADING:5 TRAILING:5 \
    SLIDINGWINDOW:4:15 MINLEN:36 &
  #java -jar $trim_dir/trimmomatic-0.35.jar SE -threads $cores -phred33 $F1 TRIM_SE/$name.trim.fq.gz ILLUMINACLIP:$trim_dir/adapters/TruSeq3-SE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36 &
done
wait

# Trim paired end reads
for F1 in `ls fastq/PE/*R1*.gz`; do
  name=${F1/.R1.gz/}; name=${name/.R1.gz}; name=$(basename $name)
#  name=${F1/_R1_001.fastq.gz/}; name=${name/.R1.fastq.gz}; name=$(basename $name)
  F2=${F1/R1/R2}
  $SRUN java -jar $trim_dir/trimmomatic-0.35.jar \
    PE -threads $cores -phred33 $F1 $F2 \
    TRIM_PE/$name.R1.trim.fq.gz TRIM_PE/$name.R1.unpaired.fq.gz \
    TRIM_PE/$name.R2.trim.fq.gz TRIM_PE/$name.R2.unpaired.fq.gz \
    ILLUMINACLIP:$trim_dir/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 \
    SLIDINGWINDOW:4:15 MINLEN:36 &
  #java -jar $trim_dir/trimmomatic-0.35.jar PE -threads $cores -phred33 $F1 $F2 TRIM_PE/$name.R1.trim.fq.gz TRIM_PE/$name.R1.unpaired.fq.gz TRIM_PE/$name.R2.trim.fq.gz TRIM_PE/$name.R2.unpaired.fq.gz ILLUMINACLIP:$trim_dir/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36 &
done
wait

echo ------------------- done TRIM_PE and TRIM_SE;echo `date +"%Y-%m-%d %T"`;echo;echo;


#----------------------------------------------------------------------------------------
# 4. Run FastQC on the trimmed reads
#----------------------------------------------------------------------------------------
if [ ! -d FASTQC_TRIM ]; then mkdir FASTQC_TRIM; fi
if [ ! -d FASTQC_TRIM/SE ]; then mkdir FASTQC_TRIM/SE; fi
if [ ! -d FASTQC_TRIM/PE ]; then mkdir FASTQC_TRIM/PE; fi

# FastQC on SE trimmed files
for f in `find $WORKING_DIR/TRIM_SE -name "*trim.fq.gz"`; do
  $SRUN fastqc --outdir FASTQC_TRIM/SE --threads $cores $f &
  #fastqc --outdir FASTQC_TRIM/SE --threads $cores $f &
done;
wait

# FastQC on PE trimmed files
for f in `find $WORKING_DIR/TRIM_PE -name "*trim.fq.gz"`; do
  $SRUN fastqc --outdir FASTQC_TRIM/PE --threads $cores $f &
  #fastqc --outdir FASTQC_TRIM/PE --threads $cores $f &
done;
wait

# Write FastQC results for each sample after Trimming to summary_FASTQC_TRIM.txt:
#    - Reads post trimming
#    - Post trimming length 
#    - %GC post trimming
cd FASTQC_TRIM/SE
for f in `ls *.zip`; do unzip $f; done
grep "Total Sequences" */fastqc_data.txt > 1_reads.txt
grep "Sequence length" */fastqc_data.txt > 1_readLength.txt
grep "^%GC" */fastqc_data.txt > 1_GC.txt
sed -i 's/ /\_/g; s/\.trim\_fastqc\/fastqc_data.txt\:/\t/g; s/\_001//g; s/\.R1//g; s/\.fastq.gz//g; s/\.fq.gz//g;'  *.txt
awk 'BEGIN{print "SampleName\tReadsPostTrimming\tPostTrimLength\t%GCPostTrim"} {a[$1]=a[$1]$3"\t"} END{for (i in a){print i "\t" a[i]}}' 1_reads.txt 1_readLength.txt 1_GC.txt > 1_temp.txt
(head -n 1 1_temp.txt && tail -n +2 1_temp.txt | sort) > summary.txt
cd ../PE
for f in `ls *.zip`; do unzip $f; done
grep "Total Sequences" *trim*/fastqc_data.txt > 1_reads.txt
grep "Sequence length" *trim*/fastqc_data.txt > 1_readLength.txt
grep "^%GC" *trim*/fastqc_data.txt > 1_GC.txt
sed -i 's/ /\_/g; s/\.trim\_fastqc\/fastqc_data.txt\:/\t/g; s/\_R1//g; s/\_001//g; s/\.R1//g; s/\.fastq//g' *.txt
awk '{sum+=$3} NR%2==1{printf $1 "\t"} NR%2==0{print sum; sum=0;}' 1_reads.txt > 1_reads2.txt
awk 'NR%2==1{print $1 "\t" $3}' 1_readLength.txt > 1_readLength2.txt
awk 'NR%2==1{print $1 "\t" $3}' 1_GC.txt > 1_GC2.txt
awk 'BEGIN{print "SampleName\tReadsPostTrimming\tPostTrimLength\t%GCPostTrim"} {a[$1]=a[$1]$2"\t"} END{for (i in a){print i "\t" a[i]}}' 1_reads2.txt 1_readLength2.txt 1_GC2.txt > 1_temp.txt
(head -n 1 1_temp.txt && tail -n +2 1_temp.txt | sort) > summary.txt
cd ..
cat SE/summary.txt > $WORKING_DIR/summary_FASTQC_TRIM.txt
tail -n +2 PE/summary.txt >> $WORKING_DIR/summary_FASTQC_TRIM.txt
rm $WORKING_DIR/FASTQC_TRIM/*/1_*.txt
rm $WORKING_DIR/FASTQC_TRIM/*/*.html
rm -r $WORKING_DIR/FASTQC_TRIM/*/*fastqc
cd ..

echo ------------------- done FASTQC_TRIM;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 5. Deduplicate (a.k.a. collapse) reads that are exact copies before mapping
#    with fastx for single-end reads and fastuniq for paired-end reads
#----------------------------------------------------------------------------------------
if [ ! -d COLLAPSED_SE ]; then mkdir COLLAPSED_SE; fi
if [ ! -d COLLAPSED_PE ]; then mkdir COLLAPSED_PE; fi

# Collapse single-end reads
#for f in `ls TRIM_SE/*.trim.fq.gz | grep -Ev "R[1|2]"`; do
for f in `ls TRIM_SE/*.trim.fq.gz`; do
  $SRUN ${SCRIPTS_DIR}/run_fastx_collapser_se.sh $f $PATH &
  #./${SCRIPTS_DIR}/run_fastx_collapser_se.sh $f $PATH &
done
wait

# Collapse paired-end reads
for F1 in `ls TRIM_PE/*.trim.fq.gz | grep R1`; do
  $SRUN ${SCRIPTS_DIR}/run_fastuniq_pe.sh $F1 $PATH &
  #./${SCRIPTS_DIR}/run_fastuniq_pe.sh $F1 $PATH &
done
wait

# Write read results for each sample after deduplicating reads to summary.txt:
#    - Trim Dedup Reads
cd COLLAPSED_SE
if [ -e summary.txt ]; then rm summary.txt; fi
ls *.trim.uniq.fq.gz | while read line; do zcat $line | echo -e "$line\t" $((`wc -l`/4)) >> summary.txt ; done
sed -i 's/\.R[1|2]//g; s/\.trim//g; s/\.uniq//g; s/\.fq//g; s/\.fastq//g; s/\.gz//g' summary.txt
awk 'BEGIN{print "SampleName\tTrimDedupReads"} NR{print $1 "\t" $2*2}' summary.txt > 1_temp.txt
(head -n 1 1_temp.txt && tail -n +2 1_temp.txt | sort) > summary.txt
cd ../COLLAPSED_PE
if [ -e summary.txt ]; then rm summary.txt; fi
#ls *.trim*.uniq*.fq* | while read line; do cat $line | echo -e "$line\t" $((`wc -l`/4)) >> summary.txt ; done
ls *.trim*.uniq*.fq*gz | while read line; do zcat $line | echo -e "$line\t" $((`wc -l`/4)) >> summary.txt ; done
sed -i 's/\.R[1|2]//g; s/\.trim//g; s/\.uniq//g; s/\.fq//g; s/\.fastq//g; s/\.gz//g' summary.txt
awk 'BEGIN{print "SampleName\tTrimDedupReads"}{sum+=$2} NR%2==1{printf $1 "\t"} NR%2==0{print sum; sum=0;}' summary.txt > 1_temp.txt
(head -n 1 1_temp.txt && tail -n +2 1_temp.txt | sort) > summary.txt
cd ..
cat COLLAPSED_SE/summary.txt > summary_Trim_Dedup.txt
tail -n +2 COLLAPSED_PE/summary.txt >> summary_Trim_Dedup.txt
rm COLLAPSED_SE/1_*.txt; rm COLLAPSED_PE/1_*.txt

echo ------------------- done COLLAPSED_SE COLLAPSED_PE;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 6. Map with BWA-MEM
#    *note: The genome must first be indexed w/ bwa index. Example for rheMac8 genome:
#          REF=$WORKING_DIR/RheMac8/rheMac8.fa
#          name=$(basename $REF); name=$WORKING_DIR/RheMac8/$name;
#          bwa index -a bwtsw $name $name.amb $name.ann $name.bwt $name.pac $name.sa
#
#    This script runs bwa-mem, and filters for quality of 30
#    Outputs are:
#      <read_name>.mapped_reads   # Read counts at each step
#      <read_name>.q30.sort.bam
#      <read_name>.q30.sort.bam.bai
#----------------------------------------------------------------------------------------
if [ ! -d MAPPED_SE ]; then mkdir MAPPED_SE; fi
if [ ! -d MAPPED_PE ]; then mkdir MAPPED_PE; fi

# Map single-end reads
#for f in `ls COLLAPSED_SE/*.fq.gz | grep -Ev "R[1|2]"`; do
for f in `ls COLLAPSED_SE/*.fq.gz`; do
  $SRUN ${SCRIPTS_DIR}/run_bwa_mem_se.sh $f $PATH $WORKING_DIR $REF &
  #./${SCRIPTS_DIR}/run_bwa_mem_se.sh $f $PATH $WORKING_DIR $REF &
done

# Map paired-end reads
for F1 in `ls COLLAPSED_PE/*.fq.gz | grep R1`;do
  $SRUN ${SCRIPTS_DIR}/run_bwa_mem_pe.sh $F1 $PATH $WORKING_DIR $REF &
  #./${SCRIPTS_DIR}/run_bwa_mem_pe.sh $F1 $PATH $WORKING_DIR $REF &
done
wait

# Write mapped results for each sample after alignment to summary_MAPPED.txt:
#    - Mapped Reads
#    - Reads Passing Quality of 30
#    - Unique Mapping Positions
$SRUN ${SCRIPTS_DIR}/summarize.sh $WORKING_DIR
#./$SCRIPTS_DIR/summarize.sh $WORKING_DIR

echo ------------------- done MAPPED_SE MAPPED_PE and stats;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 7. Convert bam files to bed files for viewing
#----------------------------------------------------------------------------------------
for f in `ls MAPPED_SE/*.q30.sort.bam`; do
  $SRUN bamToBed -i $f > ${f/.bam/.bed} &
  #bamToBed -i $f > ${f/.bam/.bed} &
done
for f in `ls MAPPED_PE/*.q30.sort.bam`; do
  $SRUN bamToBed -i $f > ${f/.bam/.bed} &
  #bamToBed -i $f > ${f/.bam/.bed} &
done
wait

echo ------------------- done bed files;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 8. Generate statistics summary table containing:
#      - sample reads
#      - trimmed reads
#      - mapped reads
#----------------------------------------------------------------------------------------
cd $WORKING_DIR

# Write number of reads mapped to chrM for all samples to chrM_counts.txt:
#   - mtDNA read counts
echo "SampleName mtDNA_ReadCounts" > chrM_counts.txt
for f in `ls $WORKING_DIR/MAPPED_*/*i.q30.sort.bam`; do
  name=${f/.q30.sort.bam/}
  name=$(basename $name)
  count=`samtools view $f | grep "$chrMT" | wc -l`
  echo $name $count >> chrM_counts.txt
done

# Write number of reads mapped to chrY for all samples to chrY_counts.txt:
#   - chrY read counts
echo "SampleName chrY_ReadCounts" > chrY_counts.txt
for f in `ls $WORKING_DIR/MAPPED_*/*.q30.sort.bam`; do
  name=${f/.q30.sort.bam/}
  name=$(basename $name)
  count=`samtools view $f | grep "$chrY" | wc -l`
  echo $name $count >> chrY_counts.txt
done

# Generate table of all summary stats before counts and copy number
(head -n 2 chrM_counts.txt && tail -n +3 chrM_counts.txt | sort -V) > chrM_counts_temp.txt
(head -n 2 chrY_counts.txt && tail -n +3 chrY_counts.txt | sort -V) > chrY_counts_temp.txt
(head -n 2 summary_FASTQC.txt && tail -n +3 summary_FASTQC.txt | sort -V) > summary_FASTQC_temp.txt
(head -n 2 summary_FASTQC_TRIM.txt && tail -n +3 summary_FASTQC_TRIM.txt | sort -V) > summary_FASTQC_TRIM_temp.txt
(head -n 2 summary_Trim_Dedup.txt && tail -n +3 summary_Trim_Dedup.txt | sort -V) > summary_Trim_Dedup_temp.txt
(head -n 2 summary_MAPPED.txt && tail -n +3 summary_MAPPED.txt | sort -V) > summary_MAPPED_temp.txt
join chrM_counts_temp.txt chrY_counts_temp.txt > 1_tempSummary1.txt
join 1_tempSummary1.txt summary_FASTQC_temp.txt > 1_tempSummary2.txt
join 1_tempSummary2.txt summary_FASTQC_TRIM_temp.txt > 1_tempSummary1.txt
join 1_tempSummary1.txt summary_Trim_Dedup_temp.txt > 1_tempSummary2.txt
join 1_tempSummary2.txt summary_MAPPED_temp.txt > 1_tempSummary1.txt 
rm summary_*_temp.txt

# Add additional statistics: 
#   - %Trimmed = (RawReads - ReadsPostTrimmed)/ReadsPostTrimming
#   - %Unique = TrimDedupRead/ReadsPostTrimming
#   - %Mapped = Mapped/TrimDedupReads
#   - %PassQ30 = PassQ30/TrimDedupReads
#   - Non-RedundantFraction = UniqueMappingPOsitions/Mapped
head -n 1 1_tempSummary1.txt > 1_tempSummaryHeader.txt
tail -n +2 1_tempSummary1.txt > 1_tempSummary2.txt
#awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, ($4-$7)/$7*100, $8, $9, $10, ($10/$7)*100, $11, $12, $13,  ($11/$10)*100, ($12/$10)*100, ($13/$11)}' 1_tempSummary2.txt > 1_tempSummary1.txt
#awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, "%Trimmed", $8, $9, $10, "%Unique", $11, $12, $13, "%Mapped", "%PassQ30", "NRF(Non-RedundantFraction)" }' 1_tempSummaryHeader.txt > 1_tempSummary2.txt
awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, ($5-$8)/$8*100, $9, $10, $11, ($11/$8)*100, $12, $13, $14, ($12/$11)*100, ($13/$11)*100, ($14/$12)}' 1_tempSummary2.txt > 1_tempSummary1.txt
awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, "%Trimmed", $9, $10, $11, "%Unique", $12, $13, $14, "%Mapped", "%PassQ30", "NRF(Non-RedundantFraction)" }' 1_tempSummaryHeader.txt > 1_tempSummary2.txt
cat 1_tempSummary2.txt > $WORKING_DIR/VNOWCHI_Summary.txt
cat 1_tempSummary1.txt >> $WORKING_DIR/VNOWCHI_Summary.txt
rm 1_tempSummary*.txt

echo ------------------- done VNOWCHI_SummaryTable;echo `date +"%Y-%m-%d %T"`;echo;


#----------------------------------------------------------------------------------------
# 9. Get counts and ratios in each bin for the 4 bin sizes: 500, 1000, 2000, 4000
#----------------------------------------------------------------------------------------
if [ ! -d COUNTS ]; then mkdir COUNTS; fi
if [ ! -d COUNTS/SE ]; then mkdir COUNTS/SE; fi
if [ ! -d COUNTS/PE ]; then mkdir COUNTS/PE; fi

# Get counts and ratios for SE
for bins in 500 1000 2000 4000; do
  if [ ! -d "$WORKING_DIR/COUNTS/SE/$bins" ]; then
     mkdir $WORKING_DIR/COUNTS/SE/$bins
  fi
  for f in `ls MAPPED_SE/*.q30.sort.bam`; do
    out=$(basename $f); out=${out/.q30.sort.bam/.counts}
    $SRUN ${SCRIPTS_DIR}/count_bin_reads.R $WORKING_DIR/FIBROBLASTS/BINS/SE/xx-euploid-1_all_merge_${bins}.bins $f $WORKING_DIR/COUNTS/SE/$bins/$out &
    #./${SCRIPTS_DIR}/count_bin_reads.R $WORKING_DIR/FIBROBLASTS/BINS/SE/xx-euploid-1_all_merge_${bins}.bins $f $WORKING_DIR/COUNTS/SE/$bins/$out &
  done
done
wait

# Get counts and ratios for PE
for bins in 500 1000 2000 4000; do
  if [ ! -d "$WORKING_DIR/COUNTS/PE/$bins" ]; then
    mkdir $WORKING_DIR/COUNTS/PE/$bins
  fi
  for f in `ls MAPPED_PE/*.q30.sort.bam`; do
    out=$(basename $f); out=${out/.q30.sort.bam/.counts}
    $SRUN ${SCRIPTS_DIR}/count_bin_reads.R $WORKING_DIR/FIBROBLASTS/BINS/PE/xx-euploid-1_all_${bins}.bins $f $WORKING_DIR/COUNTS/PE/$bins/$out &
    #./${SCRIPTS_DIR}/count_bin_reads.R $WORKING_DIR/FIBROBLASTS/BINS/PE/xx-euploid-1_all_${bins}.bins $f $WORKING_DIR/COUNTS/PE/$bins/$out &
  done
done
wait

echo ------------------- done get counts,ratios in each 4 bins;echo `date +"%Y-%m-%d %T"`;echo;



COMMENT

#----------------------------------------------------------------------------------------
# 10. Get copy number calls and aggregate data to classify embryo ploidy and sex status     
#----------------------------------------------------------------------------------------
# Get copy number calls
#   - *.cn.txt    CBS copy number and ratios used for plots
#   - *.all.txt   CBS and HMM results for all windows
#   - *.CNV.txt   Condensed CBS and HMM results with ploidy and sex status
for f in `ls $WORKING_DIR/COUNTS/*/*/*.counts`; do
  $SRUN ${SCRIPTS_DIR}/get_copy_number_and_classify.R $f $CHR_SIZE_TABLE $WORKING_DIR/VNOWCHI_Summary.txt &
  #./${SCRIPTS_DIR}/get_copy_number_and_classify.R $f $CHR_SIZE_TABLE $WORKING_DIR/VNOWCHI_Summary.txt &
done
wait

# Uncomment the script if you want the following: 
# Concatenate CBS and HMM results for all windows for all samples by SE or PE and bins
#   - all_SE_${bins}.txt  all window results for SE based on {bin}
#   - all_SE_${bins}.txt  all window results for PE based on {bin}
#for bins in 500 1000 2000 4000; do
#  touch $WORKING_DIR/COUNTS/all_SE_${bins}.txt
#  for f in `ls $WORKING_DIR/COUNTS/SE/${bins}/*.all.txt`; do
#    if [ $(wc -l < $WORKING_DIR/COUNTS/all_SE_${bins}.txt) -ge 1 ]; then
#      tail -n +2 $f >> $WORKING_DIR/COUNTS/all_SE_${bins}.txt
#    else
#      cat $f >> $WORKING_DIR/COUNTS/all_SE_${bins}.txt
#    fi
#  done
#done
#for bins in 500 1000 2000 4000; do
#  touch $WORKING_DIR/COUNTS/all_PE_${bins}.txt
#  for f in `ls $WORKING_DIR/COUNTS/PE/${bins}/*.all.txt`; do
#    if [ $(wc -l < $WORKING_DIR/COUNTS/PE/${bins}/all_${bins}.txt) -ge 1 ]; then
#      tail -n +2 $f >> $WORKING_DIR/COUNTS/all_PE_${bins}.txt
#    else
#      cat $f >> $WORKING_DIR/COUNTS/all_PE_${bins}.txt
#    fi
#  done
#done


# Concatenate condensed CBS and HMM results with ploidy and sex status for all samples by SE or PE and bins 
#   - CNV_SE_${bins}.txt  CNV results for SE based on {bin}
#   - CNV_PE_${bins}.txt  CNV results for PE based on {bin}
for bins in 500 1000 2000 4000; do
  if [ -f $WORKING_DIR/COUNTS/CNV_SE_${bins}.txt ]; then rm $WORKING_DIR/COUNTS/CNV_SE_${bins}.txt; fi
  if [ -f $WORKING_DIR/COUNTS/CNV_PE_${bins}.txt ]; then rm $WORKING_DIR/COUNTS/CNV_PE_${bins}.txt; fi
done

for bins in 500 1000 2000 4000; do
  touch $WORKING_DIR/COUNTS/CNV_SE_${bins}.txt
  for f in `ls $WORKING_DIR/COUNTS/SE/$bins/*.CNV.txt`; do
    if [ $(wc -l < $WORKING_DIR/COUNTS/CNV_SE_${bins}.txt) -ge 1 ]; then
      tail -n +2 $f >> $WORKING_DIR/COUNTS/CNV_SE_${bins}.txt
    else 
      cat $f >> $WORKING_DIR/COUNTS/CNV_SE_${bins}.txt
    fi
  done
done
for bins in 500 1000 2000 4000; do
  touch $WORKING_DIR/COUNTS/CNV_PE_${bins}.txt
  for f in `ls $WORKING_DIR/COUNTS/PE/$bins/*.CNV.txt`; do
    if [ $(wc -l < $WORKING_DIR/COUNTS/CNV_PE_${bins}.txt) -ge 1 ]; then
      tail -n +2 $f >> $WORKING_DIR/COUNTS/CNV_PE_${bins}.txt
    else
      cat $f >> $WORKING_DIR/COUNTS/CNV_PE_${bins}.txt
    fi
  done
done
wait

COMMENT

# Summarize copy number variations per embryo based on individual samples
# for SE or PE based on bins:
#  - CNV_{SE|PE}_${bins}.sampleSummary.txt  Summary of all individual samples
#  - CNV_{SE|PE}_${bins}.chrSummary.txt     Summary of whole and segmental chromosome missegregation per embryo
#  - CNV_{SE|PE}_${bins}.chrFreq.txt        Summary of whole and segmental chromosome missegregation frequency
#  - CNV_{SE|PE}_${bins}.segSummary.txt     Summary of segmental chromosome missegregation per embryo
#  - CNV_{SE|PE}_${bins}.segFreq.txt        Summary of segmental chromosome missegregation frequency
for bins in 500 1000 2000 4000; do
  f="$WORKING_DIR/COUNTS/CNV_SE_${bins}.txt"
  $SRUN ${SCRIPTS_DIR}/aggregate_cnv.R $f $GROUP_INFO_SE_TABLE &
  #./${SCRIPTS_DIR}/aggregate_cnv.R $f $GROUP_INFO_SE_TABLE &
done

for bins in 500 1000 2000 4000; do
  f="$WORKING_DIR/COUNTS/CNV_PE_${bins}.txt"
  $SRUN ${SCRIPTS_DIR}/aggregate_cnv.R $f $GROUP_INFO_PE_TABLE &
  #./${SCRIPTS_DIR}/aggregate_cnv.R $f $GROUP_INFO_PE_TABLE &
done
wait

exit


#----------------------------------------------------------------------------------------
# 11. Create plots for all samples
#----------------------------------------------------------------------------------------
# Plot CBS copy number using *.cn.txt
for f in `ls $WORKING_DIR/COUNTS/*/*/*.cn.txt`; do 
  $SRUN ${SCRIPTS_DIR}/plot_copy_number.R $f &
  #./${SCRIPTS_DIR}/plot_copy_number.R $f &
done

# Plot raw counts for the 4000 bins for all samples
for f in `ls $WORKING_DIR/COUNTS/*/4000/*.cn.txt`; do
  $SRUN ${SCRIPTS_DIR}/plot_raw.R $f &
  #./${SCRIPTS_DIR}/plot_raw.R $f &
done

# Make whole genome plots for 16 samples at a time
# Will produce mulit_<num>.png files in the given directory
$SRUN ${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/PE/500/ &
$SRUN ${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/PE/1000/ &
$SRUN ${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/PE/2000/ &
$SRUN ${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/PE/4000/ &
$SRUN ${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/SE/500/ &
$SRUN ${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/SE/1000/ &
$SRUN ${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/SE/2000/ &
$SRUN ${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/SE/4000/ &
#./${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/PE/500/ &
#./${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/PE/1000/ &
#./${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/PE/2000/ &
#./${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/PE/4000/ &
#./${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/SE/500/ &
#./${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/SE/1000/ &
#./${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/SE/2000/ &
#./${SCRIPTS_DIR}/plot_many.R $WORKING_DIR/COUNTS/SE/4000/ &
wait

echo ------------------- done copy number calls and plots made;echo `date +"%Y-%m-%d %T"`;echo;

