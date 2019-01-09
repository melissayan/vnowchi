#!/usr/bin/bash

#----------------------------------------------------------------------------
#  Download Tools for Variable Non-Overlapping Window CBS & HMM Intersect Pipeline
#----------------------------------------------------------------------------
#  Downloads tools required for Variable Non-Overlapping Window CBS & HMM 
#  Intersect Pipeline into your user ACC home directory.
# 
#  Input: your username
#	ex. $./download.sh username
#  Output: None
#
#  Software versions:
#    - fastqc: v0.10.1
#    - Trimmomatic v0.35
#    - fastx: v0.0.13
#    - BWA-mem: 0.7.9a-r786
#    - samtools: 0.1.19-44428cd
#    - Bedtools: v2.25.0
#    - fastuniq: v1.1
#
#----------------------------------------------------------------------------

USERNAME=$1
CUR_DIR=$(pwd)
USER_DIR=/home/users/$USERNAME
BIN_DIR=/home/users/$USERNAME/bin


#fastqc v0.11.5
cd $CUR_DIR
wget -e robots=off --cut-dirs=3 --reject="index.html*" --no-parent --recursive --relative --level=1 --no-directories https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip  fastqc_v0.11.5.zip -d $BIN_DIR
chmod 755 $BIN_DIR/FastQC/fastqc
ln -s $BIN_DIR/FastQC/fastqc $BIN_DIR/fastqc


#Trimmomatic v0.35
cd $CUR_DIR
wget -e robots=off --cut-dirs=3 --reject="index.html*" --no-parent --recursive --relative --level=1 --no-directories http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip
unzip Trimmomatic-0.35.zip -d $BIN_DIR
chmod 755 $BIN_DIR/Trimmomatic-0.35/trimmomatic-0.35.jar
ln -s $BIN_DIR/Trimmomatic-0.35/trimmomatic-0.35.jar $BIN_DIR/trimmomatic-0.35.jar 


#fastx v0.0.13
cd $CUR_DIR
wget -e robots=off --cut-dirs=3 --reject="index.html*" --no-parent --recursive --relative --level=1 --no-directories http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
mkdir $BIN_DIR/fastx
tar xvjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 -C $BIN_DIR/fastx
for file in $BIN_DIR/fastx/bin/fast*; do
    tool=$(echo $file | awk -F/ '{print $NF}') 
    ln -s $file $BIN_DIR/$tool
done


#BWA-mem v0.7.9a-r786
cd $CUR_DIR
wget -e robots=off --cut-dirs=3 --reject="index.html*" --no-parent --recursive --relative --level=1 --no-directories https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.9a.tar.bz2
tar xvjf bwa-0.7.9a.tar.bz2 -C $BIN_DIR
cd $BIN_DIR/bwa-0.7.9a
./configure
make
make install
ln -s $BIN_DIR/bwa-0.7.9a/bwa $BIN_DIR/bwa


#samtools: 0.1.19-44428cd
cd $CUR_DIR
wget -e robots=off --cut-dirs=3 --reject="index.html*" --no-parent --recursive --relative --level=1 --no-directories https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
tar xvjf samtools-0.1.19.tar.bz2 -C $BIN_DIR
cd $BIN_DIR/samtools-0.1.19
./configure
make
make install
ln -s $BIN_DIR/samtools-0.1.19/samtools $BIN_DIR/samtools


#Bedtools v2.25.0
cd $CUR_DIR
wget -e robots=off --cut-dirs=3 --reject="index.html*" --no-parent --recursive --relative --level=1 --no-directories https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar xvzf bedtools-2.25.0.tar.gz -C $BIN_DIR
cd $BIN_DIR/bedtools2
make
make install
ln -s $BIN_DIR/bedtools2/bin/bedtools $BIN_DIR/bedtools


#fastuniq v1.1
cd $CUR_DIR
wget -e robots=off --cut-dirs=3 --reject="index.html*" --no-parent --recursive --relative --level=1 --no-directories https://sourceforge.net/projects/fastuniq/files/FastUniq-1.1.tar.gz
tar xvzf FastUniq-1.1.tar.gz -C $BIN_DIR
cd $BIN_DIR/FastUniq/source
make
ln -s $BIN_DIR/FastUniq/source/fastuniq $BIN_DIR/fastuniq

