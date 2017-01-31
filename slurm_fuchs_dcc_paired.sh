#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -J "FUCHS"
#SBATCH -o FUCHS_shell_out.txt
#SBATCH -o FUCHS_errors.txt

# check if we have 3 arguments
if [ ! $# == 8 ]; then
  echo "Usage: $0 [Paired BAM file] [BED file] [target dir e.g. /path/to/data/] [Sample name e.g. Sample_A] [Mate 1 chimeric junction file] [Mate 2 chimeric junction file] [CircRNACount file] [Minimum read support] [Position of the exon number in BED annotation] [Minimum MAPQ] [Steps to skip]"
  exit
fi

FUCHS mock $1 $2 $3 $4 -j $5 -m $6 -c $7 --tmp /scratch/global_tmp/ -r $8 -e $9 -q $10 -sS $11 
