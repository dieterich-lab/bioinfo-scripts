#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Date:   Tuesday, May 3, 2016 6:51 PM
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @Last modified time: Wednesday, May 4, 2016 11:53 AM
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=1G
#SBATCH -J "FastQC"
#SBATCH -c 1
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de


# we are good with 1GB, the manual says it allocates ~256MB per thread

module load fastqc

# check if we have 3 arguments
if [ ! $# == 2 ]; then
  echo $#
  echo "Usage: $0 [Read file] [output directory]"
  exit
fi

#fastqc $1 -o $2 -a /biosw/flexbar/Adapter.fa -d /scratch/
fastqc $1 -o $2 -d /scratch/global_tmp/
