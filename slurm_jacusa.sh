#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 40
#SBATCH --mem=250G
#SBATCH -J "jacusa"
#SBATCH -p long
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=tobias.jakobi@med.uni-heidelberg.de

module load java

export set1=$1
export set2=$2
export output=$3

#java -Xmx250g -jar /biosw/jacusa/1.2.3/build/JACUSA_v1.2.3.jar call-2 -s -c 5 -P RF-FIRSTSTRAND,RF-FIRSTSTRAND -p 20 -W 1000000 --filterNM_1 5 --filterNM_2 5 -F 1024 -a D,M,Y -r ${output} ${set1} ${set2}

java -Xmx250g -jar /biosw/jacusa/1.2.4/JACUSA_v1.2.4.jar call-2 -s -c 5 -P RF-FIRSTSTRAND,RF-FIRSTSTRAND -p 10 -W 1000000 -F 1024 --filterNM_1 5 --filterNM_2 5 -T 1 -a D,M,Y -u DirMult:showAlpha -R -r ${output} ${set1} ${set2}
￼
