#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Date:   Tuesday, May 3, 2016 6:51 PM
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @Last modified time: Wednesday, May 4, 2016 10:36 AM
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=4G
#SBATCH -J flexbar paired


# check if we have 3 arguments
if [ ! $# == 3 ]; then
  echo "Usage: $0 [Read 1 file] [Read 2 file] [target dir e.g. /tmp/]"
  exit
fi

# $1 -> Read 1
# $2 -> Read 3
# $3 -> Target directory

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${1/_R1/} : '\(.*\)\..*\.'`target=`expr ${1/_R1/} : '\(.*\)\..*\.'`

# load the flexbar module
module load flexbar

# run on 20 CPUs
# compress with bz2
# only 30nt or longer
# no uncalled bases
# quality min phred 28
# use sanger quality values (i.e. Illumina 1.9+ encoding)
flexbar -r $1 -p $2 -t $3/$target  -n 20 -z BZ2 -m 30 -u 0 -q 28 -a /biosw/flexbar/Adapter.fa -f sanger
