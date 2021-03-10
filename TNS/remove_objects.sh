#!/usr/bin/env bash

grep 'Update' reclass_list| awk -F: '{print $1}' | awk '{print $3}' | while read AT; do
    for i in $(ls */*txt); do
	echo $AT $i
	grep -v ${AT} ${i} > tmp
	diff tmp $i
	mv tmp $i
    done
done
