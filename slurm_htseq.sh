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
#SBATCH -c 1
#SBATCH --mem=10G
#SBATCH -J "htseq"


# check if we have 3 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [Read 1 file] [gff] [target dir e.g. /tmp/]"
  exit
fi

# $1 -> Genome index
# $2 -> Read 1
# $3 -> Read 2
# $4 -> Target directory

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${2/$5/} : '\(.*\)\..*\.'`

 htseq-count -f bam $4/$target/Aligned.sortedByCoord.out.bam /biodb/genomes/mus_musculus/GRCm38_85/GRCm38.85.gtf > $4/$target/${target}.counts
 htseq-count stranded=no -f bam $4/$target/Aligned.sortedByCoord.out.bam /biodb/genomes/mus_musculus/GRCm38_85/GRCm38.85.gtf > $4/$target/${target}.counts.unstranded
