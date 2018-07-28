#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40
#SBATCH --mem=250G
#SBATCH -J "DCC"

# check if we have 3 arguments
if [ ! $# >= 8 ]; then
  echo "Usage: $0 [Sample sheet file] [GTF file] [Genome FASTA file] [BAM list file] [Repeats] [target dir e.g. project/] [flags]"
  exit
fi

DCC @$1 $7 -D -an $2 -A $3 -B @$4 -Nr 2 2 -fg -G -k -O $6 --temp ${6}_DCC_temp/ -F -L 20 -T 40 -ss
