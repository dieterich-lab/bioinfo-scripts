#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40
#SBATCH --mem=250G
#SBATCH -J "circtools detect"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 3 arguments
if [ ! $# == 6 ]; then
  echo "Usage: $0 [Sample sheet file] [GTF file] [Genome FASTA file] [BAM list file] [Repeats] [target dir e.g. project/]"
  exit
fi

/usr/bin/time -v circtools detect @$1 -D -an $2 -A $3 -B @$4 -R $5 -M -Nr 2 2 -fg -G -k -O $6 -t ${6}_DCC_temp/ -F -L 20 -T 40
