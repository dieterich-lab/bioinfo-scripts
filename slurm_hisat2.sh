#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Date:   Friday, May 6, 2016 4:10 PM
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @Last modified time: Friday, May 6, 2016 4:22 PM
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40
#SBATCH --mem=250G
#SBATCH -J "HISAT2 alignment"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module load hisat2

# check if we have 6 arguments
if [ ! $# == 6 ]; then
  echo "Usage: $0 [HISAT2 index] [Read 1 file] [Read 2 file] [target dir e.g. /tmp/] [Read 1 marker, e.g. R1] [GTF file]"
  exit
fi

# $1 -> Genome index
# $2 -> Read 1
# $3 -> Read 2
# $4 -> Target directory

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${2/$5/} : '\(.*\)\..*\.'`

# create the target directory, HISAT2 will not do that for us
mkdir -pv $4/$target

OLD_PATH=`pwd`

hisat2 --time -p 40 -x ${1} -1 ${2} -2 ${3} -S $4/$target/Aligned.out.sam \
       --summary-file \
       --new-summary \
       --un-conc-gz $4/$target/ \
       --novel-splicesite-outfile $4/$target/novel_splice_sites.csv \
       --known-splicesite-infile ${6} \
       --rna-strandness FR

cd  $4/$target/

samtools view -bS Aligned.out.sam | samtools sort -@ 20 -m 10G -T tempo -o Aligned.out.bam /dev/stdin 
samtools index Aligned.out.bam

rm -f Aligned.out.sam

