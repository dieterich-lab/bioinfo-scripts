#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=250G
#SBATCH -J "picard-tools"


# check if we have 2 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [Input BAM] [Output BAM]"
  exit
fi

#module load picard-tools

java -Xmx95g -jar /biosw/picard-tools/2.5.0/picard.jar MarkDuplicates I=$1 O=$2 M=${1}.duplicate.metrics ASSUME_SORT_ORDER=coordinate

