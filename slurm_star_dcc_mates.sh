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
#SBATCH -c 5 
#SBATCH --mem=200G
#SBATCH -J "STAR genome alignment"

module load star

# check if we have 5 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [STAR index] [Read file] [target dir e.g. /tmp/] [GTF file] [Mate ID e.g. 1]"
  exit
fi

# $1 -> Genome index
# $2 -> Read 1
# $3 -> Read 2

target=`expr ${2/$5/} : '\(.*\)\..*\.'`

# create the target directory, STAR will not do that for us

folder=$3/${target}_mate$5


echo $folder
mkdir -pv $folder

STAR   --readFilesCommand zcat\
       --runThreadN 5\
       --genomeDir $1\
       --outSAMtype BAM Unsorted SortedByCoordinate\
       --readFilesIn $2\
       --outFileNamePrefix $folder/\
       --outReadsUnmapped Fastx \
       --quantMode GeneCounts\
       --genomeLoad NoSharedMemory\
       --outWigType bedGraph\
       --outReadsUnmapped Fastx\
       --outSJfilterOverhangMin 15 15 15 15\
       --alignSJoverhangMin 15\
       --alignSJDBoverhangMin 10\
       --outFilterMultimapNmax 20\
       --outFilterScoreMin 1\
       --outFilterMismatchNmax 999\
       --outFilterMismatchNoverLmax 0.05\
       --outFilterMatchNminOverLread 0.7\
       --alignIntronMin 20\
       --alignIntronMax 1000000\
       --alignMatesGapMax 1000000\
       --chimSegmentMin 15\
       --chimScoreMin 15\
       --chimScoreSeparation 10\
       --chimJunctionOverhangMin 15\
       --twopassMode Basic\
       --alignSoftClipAtReferenceEnds No\
       --outSAMattributes NH HI AS nM NM MD jM jI XS\
       --sjdbGTFfile $4
