#!/bin/bash

gff=$1


grep -v 'softclipped_support=0' $gff | perl -pe 's/(.+softclipped_pos=\()(\d+), (\d+)(\).+)/\2	\3	\1\2,\3\4/' | awk 'BEGIN{FS="\t";OFS="\t"};{print $3,$4,$5,$1,$2,$8,$9,$10,$11}' > ${gff/.gff/}.SC_COORDS.gff
