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
#SBATCH -c 20
#SBATCH --mem=200G
#SBATCH -J "circtools alignment"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module load star

# check if we have 6 arguments
if [ ! $# == 6 ]; then
  echo "Usage: $0 [STAR index] [Read 1 file] [Read 2 file] [target dir e.g. /tmp/] [Read 1 marker, e.g. R1] [GTF file]"
  exit
fi

# $1 -> Genome index
# $2 -> Read 1
# $3 -> Read 2
# $4 -> Target directory

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${2/$5/} : '\(.*\)\..*\.'`

# create the target directory, STAR will not do that for us
mkdir -pv $4/$target
mkdir -pv $4/$target/mate1/
mkdir -pv $4/$target/mate2/

mkdir -pv $4/$target/mate2_new/
mkdir -pv $4/$target/mate1_new/

STAR --runThreadN 20 \
       --genomeDir $1 \
       --outSAMtype BAM SortedByCoordinate \
       --readFilesIn $2 $3 \
       --readFilesCommand zcat \
       --outFileNamePrefix $4/$target/ \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15 \
       --chimOutType Junctions SeparateSAMold\
       --sjdbGTFfile $6 

gzip $4/$target/Unmapped.out.mate1
gzip $4/$target/Unmapped.out.mate2

STAR --runThreadN 20 \
       --genomeDir $1 \
       --outSAMtype BAM SortedByCoordinate \
       --readFilesIn $4/$target/Unmapped.out.mate1.gz \
       --readFilesCommand zcat \
       --outFileNamePrefix $4/$target/mate1/ \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --seedSearchStartLmax 30 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15 \
       --chimOutType Junctions SeparateSAMold\
       --sjdbGTFfile $6


STAR --runThreadN 20 \
       --genomeDir $1 \
       --outSAMtype BAM SortedByCoordinate \
       --readFilesIn $4/$target/Unmapped.out.mate2.gz \
       --readFilesCommand zcat \
       --outFileNamePrefix $4/$target/mate2/ \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --seedSearchStartLmax 30 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15 \
       --chimOutType Junctions SeparateSAMold\
       --sjdbGTFfile $6

## map from scratch

STAR --runThreadN 20 \
       --genomeDir $1 \
       --outSAMtype BAM SortedByCoordinate \
       --readFilesIn $2 \
       --readFilesCommand zcat \
       --outFileNamePrefix $4/$target/mate1_new/ \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --seedSearchStartLmax 30 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15 \
       --chimOutType Junctions SeparateSAMold\       
       --sjdbGTFfile $6


STAR --runThreadN 20 \
       --genomeDir $1 \
       --outSAMtype BAM SortedByCoordinate \
       --readFilesIn $3 \
       --readFilesCommand zcat \
       --outFileNamePrefix $4/$target/mate2_new/ \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --seedSearchStartLmax 30 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15 \
       --chimOutType Junctions SeparateSAMold\
       --sjdbGTFfile $6
