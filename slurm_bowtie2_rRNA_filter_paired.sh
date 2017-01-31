#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Date:   Wednesday, May 4, 2016 11:14 AM
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @Last modified time: Friday, May 6, 2016 4:17 PM
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -J "bowtie2 rRNA filtering"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 5 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [rRNA index argument] [Read 1 file] [Read 2 file] [target dir e.g. /awesome/project/] [R1 marker, e.g. R1 or 1_sequence]"
  exit
fi

# $1 -> rRNA index
# $2 -> Read 1
# $3 -> Read 2
# $4 -> Target directory
# $5 -> R1 marker

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${2/$5/} : '\(.*\)\..*\.'`

# load the bowtie2 module
module load bowtie2

# SAM output goes to /dev/null
# run on 20 CPUs
# set fixed seed
# memory mapped IO for multiple instances
# display timing information
# write gz unmapping reads [== no rRNA] to target dir

bowtie2 -x $1 -1 $2 -2 $3 -S $4/$target.sam --no-unal --omit-sec-seq --threads 20 --mm --seed 1337 --time --un-conc-gz $4/$target.fastq.gz 2> $4/$target.log
