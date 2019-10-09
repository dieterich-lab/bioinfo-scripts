#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH -p small
#SBATCH --mem=50G
#SBATCH -J "circtools enrich"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 8 arguments
if [ ! $# == 9 ]; then
  echo "Usage: $0 [Chromosome sizes file] [Genome GTF file] [CLIP peak BED file] [CircRNA BED files] [output prefix] [output path] [# iterations] [TMP dir] [White list file]"
  exit
fi

#/usr/bin/time -v python3 /home/tjakobi/repos/dieterichlab/circtools/circtools/circtools.py enrich -c ${4} -b ${3} -a ${2} -g ${1} -i ${7} -I exon -p 10 -P 1 -T 1 -o ${6}/ -F ${5} -t ${8}/ -W ${9}


/usr/bin/time -v python3 /home/tjakobi/repos/dieterichlab/circtools/circtools/circtools.py enrich -c ${4} -b ${3} -a ${2} -g ${1} -i ${7} -p 10 -P 1 -T 1 -o ${6}/ -F ${5} -t ${8}/ -W ${9}

#/usr/bin/time -v circtools enrich -c ${4} -b ${3} -a ${2} -g ${1} -i ${7} -I exon -p 20 -P 1 -T 1 -o ${6}/ -F ${5} -t ${8}/
