#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH -J "prepare-dcc"
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

# check if we have 2 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [STAR source dir] [DCC destination dir]"
  exit
fi


SRC=${1}
DEST=${2}

if [ ! -d "$SRC" ]; then
  echo "Source directory $SRC does not exist!"
  exit;
fi

if [ ! -d "$DEST" ]; then
  echo "DCC directory $DEST does not exist!"
  exit;
fi

cd $SRC/

parallel ln -s `pwd`/{1}/Chimeric.out.junction ../$DEST/{1}.Chimeric.out.junction ::: * 
parallel ln -s `pwd`/{1}/Aligned.noS.bam ../$DEST/{1}.bam ::: * 
parallel ln -s `pwd`/{1}/Aligned.noS.bam.bai ../$DEST/{1}.bam.bai ::: *
parallel ln -s `pwd`/{1}/SJ.out.tab ../$DEST/{1}.SJ.out.tab ::: *

cd ../$DEST/

ls | grep bam | grep -v bai | grep -v mate | grep -v bam_files > bam_files.txt
ls | grep Chimeric.out.junction | grep -v mate  > samplesheet

