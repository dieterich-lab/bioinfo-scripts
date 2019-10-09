#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2 
#SBATCH --mem=1G
#SBATCH -J "prepare-subread"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

if [ ! $# == 2 ]; then
  echo "Usage: $0 [STAR folder] [BAM folder] "
  exit
fi

mkdir -pv ${2}
cd ${1}

parallel --verbose -j1 ln -v -s ../${1}/{1}/Aligned.noS.bam ../${2}/{1}.bam ::: *
