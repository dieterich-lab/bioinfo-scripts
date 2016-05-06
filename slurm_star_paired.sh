#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Date:   Friday, May 6, 2016 4:10 PM
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @Last modified by:   tjakobi
# @Last modified time: Friday, May 6, 2016 4:15 PM
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



STARlong --genomeDir /data/projects/departments/Christoph_Dieterich/circRNA_evolution/tobias/genomes/{=s/_.*//=}/\
--runThreadN 40\
--readFilesIn /data/projects/departments/Christoph_Dieterich/circRNA_evolution/tobias/names/{..}.1.gz\
/data/projects/departments/Christoph_Dieterich/circRNA_evolution/tobias/names/{..}.2.gz\
--sjdbGTFfile /data/projects/departments/Christoph_Dieterich/circRNA_evolution/tobias/genomes/{=s/_.*//=}/{=s/_.*//=}.gtf\
--readFilesCommand zcat --outFileNamePrefix /data/projects/departments/Christoph_Dieterich/circRNA_evolution/tobias/STAR/{...}/\
--outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \& ::: *.1.gz









# load the bowtie2 module
module load bowtie2

# SAM output goes to /dev/null
# run on 20 CPUs
# set fixed seed
# memory mapped IO for multiple instances
# display timing information
# write bz2 unmapping reads [== no rRNA] to target dir

bowtie2 -x $1 -1 $2 -2 $3 -S /dev/null --threads 20 --mm --seed 1337 --time --un-conc-bz2 $4
