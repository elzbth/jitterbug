#!/bin/bash

unfiltered_dir=$1
filtered_dir=$2
prefix=$3

cat $filtered_dir/* > $prefix.all_filtered.gff3

mergeBed -i $prefix.all_filtered.gff3 >  $prefix.all_filtered.merged.gff3

for file in $unfiltered_dir/*.gff3
do name=$(basename $file)
name=${name/.gff3/}
intersectBed -a $file -b $prefix.all_filtered.merged.gff3  -wa > $prefix.$name.recovered.gff3
echo "*"
done


