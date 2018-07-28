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
# mkdir -pv $4/$target/mate1/
# mkdir -pv $4/$target/mate2/


STAR --readFilesCommand zcat\
	--runThreadN 5\
       --genomeDir $1\
       --readFilesIn $2 $3\
       --outFileNamePrefix  $4/$target/\
       --quantMode GeneCounts\
       --genomeLoad NoSharedMemory\
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
       --sjdbGTFfile $6

# gzip $4/$target/Unmapped.out.mate1
# gzip $4/$target/Unmapped.out.mate2
# 
# 
#   STAR --readFilesCommand zcat\
#           --runThreadN 5\
#          --genomeDir $1\
#          --outSAMtype BAM SortedByCoordinate\
#          --readFilesIn $4/$target/Unmapped.out.mate1.gz \
#          --outFileNamePrefix $4/$target/mate1/\
#          --quantMode GeneCounts\
#          --genomeLoad NoSharedMemory\
#          --outReadsUnmapped Fastx\
#          --outSJfilterOverhangMin 15 15 15 15\
#          --alignSJoverhangMin 15\
#          --alignSJDBoverhangMin 10\
#          --outFilterMultimapNmax 20\
#          --outFilterScoreMin 1\
#          --outFilterMismatchNmax 999\
#          --outFilterMismatchNoverLmax 0.05\
#          --outFilterMatchNminOverLread 0.7\
#          --alignIntronMin 20\
#          --alignIntronMax 1000000\
#          --alignMatesGapMax 1000000\
#          --chimSegmentMin 15\
#          --chimScoreMin 15\
#          --chimScoreSeparation 10\
#          --chimJunctionOverhangMin 15\
#          --twopassMode Basic\
#          --alignSoftClipAtReferenceEnds No\
#          --outSAMattributes NH HI AS nM NM MD jM jI XS\
#          --sjdbGTFfile $6
#   
# 
# STAR --readFilesCommand zcat\
#           --runThreadN 5\
#          --genomeDir $1\
#          --outSAMtype BAM SortedByCoordinate\
#          --readFilesIn $4/$target/Unmapped.out.mate2.gz \
#          --outFileNamePrefix $4/$target/mate2/\
#          --quantMode GeneCounts\
#          --genomeLoad NoSharedMemory\
#          --outReadsUnmapped Fastx\
#          --outSJfilterOverhangMin 15 15 15 15\
#          --alignSJoverhangMin 15\
#          --alignSJDBoverhangMin 10\
#          --outFilterMultimapNmax 20\
#          --outFilterScoreMin 1\
#          --outFilterMismatchNmax 999\
#          --outFilterMismatchNoverLmax 0.05\
#          --outFilterMatchNminOverLread 0.7\
#          --alignIntronMin 20\
#          --alignIntronMax 1000000\
#          --alignMatesGapMax 1000000\
#          --chimSegmentMin 15\
#          --chimScoreMin 15\
#          --chimScoreSeparation 10\
#          --chimJunctionOverhangMin 15\
#          --twopassMode Basic\
#          --alignSoftClipAtReferenceEnds No\
#          --outSAMattributes NH HI AS nM NM MD jM jI XS\
#          --sjdbGTFfile $6
