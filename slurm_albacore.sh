#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40 
#SBATCH --mem=250G
#SBATCH -J "albacore"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de
#SBATCH -p long

#module unload python3*
module load albacore

# check if we have 2 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [Input folder (recursive)] [Output folder] [Kit] "
  exit
fi

# we set the flow cell to 106:
#Christoph Dieterich @cdieterich 5:14 PM
#Aha, in that case you can practically always start with FLO-MIN106... 107 is for that 1D squared chemistry which we never use

read_fast5_basecaller.py -k ${3} -f FLO-MIN106 -s ${2} -o fastq -r -q 0 -t 40 -i ${1}
