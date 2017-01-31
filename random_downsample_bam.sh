#!/bin/bash

# @Author: Tobias Jakobi <tjakobi>
# @Email:  tobias.jakobi@med.uni-heidelberg.de
# @Project: University Hospital Heidelberg, Section of Bioinformatics and Systems Cardiology
# @License: CC BY-NC-SA

# check if we have 3 arguments
if [ ! $# == 2 ]; then
  echo "Usage: $0 [BAM file] [Maximal read number]"
  exit
fi

INPUT=$1
MAX=$2

# get current read number

CURRENT=`samtools flagstat $INPUT | head -1 | cut -f 1 -d ' '`

echo $INPUT has $CURRENT reads

if [ "$CURRENT" -gt "$MAX" ]; then

	RATIO=`echo print $MAX / $CURRENT. | python`
	RATIO=${RATIO:0:4}
	echo The ratio is $RATIO
	samtools view -s $RATIO -b -o ${INPUT}.new $INPUT
else 	
	echo Library already below given read/pait threshold, creating symlink insteadt.
	ln -s ${INPUT} ${INPUT}.new
fi
