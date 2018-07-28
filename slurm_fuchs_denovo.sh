#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40
#SBATCH --mem=100G
#SBATCH -J "FUCHS de-novo"
#SBATCH -p general
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 3 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [Input directory] [Sample name]"
  exit
fi

#guided_denovo_circle_structure_parallel -c 40 -A ${1} -I ${2} -N ${3} 
#guided_denovo_circle_structure_parallel -A ${1} -c 1 -I ${2} -N ${3} -T /scratch/global_tmp/ 
#/home/tjakobi/tmp/guided_denovo_circle_structure_parallel.py -c 40 -I ${1} -N ${2} -T /scratch/global_tmp/ 
guided_denovo_circle_structure_parallel -c 40 -I ${1} -N ${2} -T /scratch/global_tmp/ 
