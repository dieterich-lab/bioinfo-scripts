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
#SBATCH -c 40
#SBATCH --mem=200G
#SBATCH -J "STAR genome alignment"


# check if we have 3 arguments
if [ ! $# == 4 ]; then
  echo "Usage: $0 [STAR index] [Read 1 file] [Read 2 file] [target dir e.g. /tmp/]"
  exit
fi

# $1 -> Genome index
# $2 -> Read 1
# $3 -> Read 2
# $4 -> Target directory

# remove the file extension and potential "R1" markings
# (works for double extension, e.g. .fastq.gz)
target=`expr ${2/_R1/} : '\(.*\)\..*\.'`



STARlong --genomeDir \
--runThreadN 40\ # 40CPUs
--readFilesIn $1 $2\
--sjdbGTFfile mes/{=s/_.*//=}/{=s/_.*//=}.gtf\
--readFilesCommand bzcat\ # use bz2 input files
--outFileNamePrefix $4\
--outSAMtype BAM SortedByCoordinate
--outSAMstrandField intronMotif
