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
#SBATCH --mem=200G
#SBATCH -J "STAR genome alignment"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

#module load star
#module load subread

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
mkdir $4/$target
STAR --genomeDir $1 --runThreadN 40 --readFilesIn $2 $3 --sjdbGTFfile $6 --readFilesCommand zcat --outFileNamePrefix $4/$target/ --outSAMtype BAM Unsorted SortedByCoordinate --outSAMstrandField intronMotif 

# --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0

#STAR   --readFilesCommand zcat\
#       --runThreadN 20\
#       --genomeDir $1\
#       --outSAMtype BAM Unsorted SortedByCoordinate\
#       --readFilesIn $2 $3\
#       --outFileNamePrefix $4/$target/\
#       --quantMode GeneCounts\
#       --genomeLoad NoSharedMemory\
#       --outWigType bedGraph\
#       --outReadsUnmapped Fastx\
#       --outSJfilterOverhangMin 15 15 15 15\
#       --alignSJoverhangMin 15\
#       --alignSJDBoverhangMin 10\
#       --outFilterMultimapNmax 20\
#       --outFilterScoreMin 1\
#       --outFilterMismatchNmax 999\
#       --outFilterMismatchNoverLmax 0.05\
#       --outFilterMatchNminOverLread 0.7\
#       --alignIntronMin 20\
#       --alignIntronMax 1000000\
#       --alignMatesGapMax 1000000\
#       --chimSegmentMin 15\
#       --chimScoreMin 15\
#       --chimScoreSeparation 10\
#       --chimJunctionOverhangMin 15\
#       --twopassMode Basic\
#       --alignSoftClipAtReferenceEnds No\
#       --outSAMattributes NH HI AS nM NM MD jM jI XS\
#       --sjdbGTFfile $6

# count features with subread package
# featureCounts -f -p -T10 -s2 -O -o ${4}/${target}/${target}.featureCounts -a $6 ${4}/${target}/Aligned.sortedByCoord.out.bam

# directly run htseq count 
#htseq-count -f bam $4/$target/Aligned.sortedByCoord.out.bam /biodb/genomes/mus_musculus/GRCm38_85/GRCm38.85.gtf > $4/$target/${target}.counts
#htseq-count stranded=no -f bam $4/$target/Aligned.sortedByCoord.out.bam /biodb/genomes/mus_musculus/GRCm38_85/GRCm38.85.gtf > $4/$target/${target}.counts.unstranded
