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
#SBATCH -c 2
#SBATCH --mem=40G
#SBATCH -J "dexseq-count"

# check if we have 3 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [Annotation] [BAM] [target dir e.g. /tmp/]"
  exit
fi

# $1 -> Annotation
# $2 -> input BAM
# $3 -> Target directory
# $4 -> yes/no/reverse strand

#param=""

if [ "$4" == "reverse" ]; then
  param="rev"
fi

python /home/tjakobi/.R/3.5/DEXSeq/python_scripts/dexseq_count.py -p yes -s ${4} -f bam -r pos  $1 $2 $3/$2.$param.txt


#python /home/tjakobi/.R/3.5/DEXSeq/python_scripts/dexseq_count.py -p yes -s reverse -f bam -r pos  $1 $2 $3/$2.rev.txt
