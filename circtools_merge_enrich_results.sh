#!/bin/bash

# Copyright (C) 2018 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# check if we have 2 arguments
if [ ! $# == 5 ]; then
  echo "Usage: $0 [Input path to enrich results] [Number of iterations (for regex)] [target dir e.g. /tmp/] [Text for plotting] [File prefix]"
  exit
fi

ITER=$2

cd $1/

awk -F '\t' '{if ($6>0) {print FILENAME"\t"$0}}' *.csv | sed 's/_.*500.*.csv//g' | grep -v circRNA_host_gene > $3/$5.csv

/home//tjakobi/repos/dieterichlab/circtools/scripts/circtools_enrich_visualization.R $3/$5.csv 0.05 10 10 $3/$5.pdf "$4" colour False
