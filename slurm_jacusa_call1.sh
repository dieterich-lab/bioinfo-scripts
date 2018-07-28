#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @License: CC BY-NC-SA

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=110G
#SBATCH -J "sbatch"
#SBATCH -p long

module load java

export set1=$1
export output=$2

echo java -Xmx100g -jar /biosw/jacusa/1.2.0/build/JACUSA_v1.2.0.jar call-1 -s -c 5 -P FR-SECONDSTRAND -p 10 -W 1000000 -F 1024 --filterNM_1 5 -T 1 -a D,M,Y -r ${output} ${set1}
java -Xmx100g -jar /biosw/jacusa/1.2.0/build/JACUSA_v1.2.0.jar call-1 -c 5 -P FR-SECONDSTRAND -p 10 -W 1000000 -F 1024 --filterNM_1 5 -T 1 -a D,M,Y -r ${output} ${set1}