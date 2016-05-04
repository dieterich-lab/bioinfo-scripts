#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Date:   Wednesday, May 4, 2016 11:14 AM
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @Last modified time: Wednesday, May 4, 2016 11:56 AM
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=4G
#SBATCH -J "bowtie2 rRNA filtering"


# check if we have 3 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [rRNA index argument] [Read 1 file] [Read 2 file] [target dir e.g. /tmp/]"
  exit
fi

# $1 -> Read 1
# $2 -> Read 3
# $3 -> Target directory

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${1/_R1/} : '\(.*\)\..*\.'`



# load the bowtie2 module
module load bowtie2

# SAM output goes to /dev/null
# run on 20 CPUs
# set fixed seed
# memory mapped IO for multiple instances
# display timing information
# write bz2 unmapping reads [== no rRNA] to target dir

bowtie2 -x $1 -1 $2 -2 $3 -S /dev/null --threads 20 --mm --seed 1337 --time --un-conc-bz2 $3
