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
#SBATCH -o DCC_shell_out.txt
#SBATCH -o DCC_errors.txt

# check if we have 3 arguments
if [ ! $# == 8 ]; then
  echo "Usage: $0 [Sample sheet file] [GTF file] [Genome FASTA file] [Mate 1 file] [Mate 2 file] [BAM list file] [Repeats] [target dir e.g. project/]"
  exit
fi

DCC @$1 -ss -D -an $2 -A $3 -Pi -mt1 @$4 -mt2 $5 -B @$6 -R $7 -M -Nr 2 2 -fg -G -k -O $8 --temp ${8}_DCC_temp/ -F -L 20 -T 40
