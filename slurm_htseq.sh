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

# $1 -> Annotation
# $2 -> input BAM
# $3 -> Target directory

htseq-count stranded=yes -f bam $2 $1 > $3/$2.counts

#htseq-count stranded=yes -f bam $4/$target/Aligned.sortedByCoord.out.bam /biodb/genomes/mus_musculus/GRCm38_85/GRCm38.85.gtf > $4/$target/${target}.counts.unstranded
