#!/bin/bash

"""
This script summarize repeat position along chromosomes.
@author: Xiaolu Wei (xiaolu_wei@urmc.rochester.edu)
"""

workingdir="TEs"
cd $workingdir || ( echo "can't open the $workingdir"; exit )
echo "Working in dir: $workingdir"

chromosome_size="dmel_scaffold2_plus0310.sizes"

targets=('G2'
		'Jockey-3'
		'DOC2'
		'G6'
		'G'
		'R1'
		'PROTOP'
		'HETA'
		'TAHRE'
		'TART'
		'BARI'
		'Jockey-1'
		)


#pull out repeats in each contig
while chromosome_size='' read -r "reads" || [[ -n "$reads" ]]; do

	contig=$(echo "$reads" | awk -F'\t' '{print $1}')	
	output=$workingdir/summary/"${contig}.summary"
	echo "working on: $output"
	
	for target in "${targets[@]}"
	do 
		workingfile=$workingdir/$target/"${target}_plot.summary"
		grep -w "$contig" "$workingfile" >> "$output"
		echo "summarizing $target"	
	done
	
	#check size of file, delete empty files
	s=$(wc -c < "$output")
	if (( $s == "0" )); then	
		rm -f "$output"
		echo "empty file, deleted"
	else
		echo "move on..."
	fi
		
done < "$chromosome_size"
echo "pulled out repeats in each contig."


#pull out repeats in each chromosome. 
#Need chromosome.txt file for each chromosome, which contains all the contig names that belong to this chromosome
chromosomes=('chromosome2'
		'chromosome3'
		'chromosome4'
		'chromosomeX'
		'chromosomeY'
		)

for chromosome in "${chromosomes[@]}"
do
	chromosome_file=$workingdir/"${chromosome}.txt"
	output_file=$workingdir/summary/"${chromosome}.summary"
	echo "working on $chromosome"
	 
	while chromosome_file='' read -r "reads" || [[ -n "$reads" ]]; do
		cat $workingdir/summary/"${reads}.summary" >> "$output_file"
		echo "cating $reads now"		
	done < "$chromosome_file"
done

echo "DONE!"

