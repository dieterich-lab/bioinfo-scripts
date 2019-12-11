#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -J "multiqc"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de


# check if we have 6 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [STAR] [FASTQC] [FLEXBAR] [BOWTIE]"
  exit
fi


multiqc $1 $2 $3 $4 -f
