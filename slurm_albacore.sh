#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 60 
#SBATCH --mem=950G
#SBATCH -J "albacore"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

#module unload python3*
module load albacore

# check if we have 2 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [input] [output]"
  exit
fi

read_fast5_basecaller.py -k SQK-RNA001 -f FLO-MIN107 -s ${2} -o fastq,fast5 -r -t 60 -i ${1}
