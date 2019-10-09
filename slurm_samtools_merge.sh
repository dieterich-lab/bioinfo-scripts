#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=50G
#SBATCH -J "samtools merge"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 6 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [output file] [input BAM suffix, e.g. bla.bam -> *bla.bam]"
  exit
fi

time samtools merge -@40 ${1} *${2}
