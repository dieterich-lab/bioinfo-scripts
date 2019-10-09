#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 15
#SBATCH --mem=20G
#SBATCH -J "flexbar single"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 3 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [Read 1 file] [target dir e.g. /tmp/] [R1 marker, e.g. _R1]"
  exit
fi

# $1 -> Read 1
# $2 -> Target directory


# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=$3

# load the flexbar module
module load flexbar

# run on 10 CPUs
# compress with bz2
# only 30nt or longer
# no uncalled bases
# quality min phred 28
# use sanger quality values (i.e. Illumina 1.9+ encoding)
flexbar -r $1 -t $2/$target  -n 15 -z GZ -m 30 -u 0  -q TAIL -qt 28 -as "AGATCGGAAGAG" -qf sanger -j
