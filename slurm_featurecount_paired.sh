#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=50G
#SBATCH -J "featureCount"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module load star
module load subread

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

# count features with subread package

# normal mode, no multi mapping reads
# featureCounts -f -p -T10 -s2 -O -o ${4}/${target}/${target}.featureCounts -a $6 ${4}/${target}/Aligned.sortedByCoord.out.bam

# run again but with multimapping reads
featureCounts -M -f -p -T10 -s1 -O -o ${4}/${target}/${target}.multi1.featureCounts -a $6 ${4}/${target}/Aligned.sortedByCoord.out.bam

# directly run htseq count 
#htseq-count -f bam $4/$target/Aligned.sortedByCoord.out.bam /biodb/genomes/mus_musculus/GRCm38_85/GRCm38.85.gtf > $4/$target/${target}.counts
#htseq-count stranded=no -f bam $4/$target/Aligned.sortedByCoord.out.bam /biodb/genomes/mus_musculus/GRCm38_85/GRCm38.85.gtf > $4/$target/${target}.counts.unstranded
