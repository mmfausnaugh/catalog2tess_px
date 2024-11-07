#!/usr/bin/env bash

#update on Nov 7, 2024, so that it only remakes the file if the transient is in there
grep 'Update' reclass_list| awk -F: '{print $1}' | awk '{print $3}' | while read AT; do
    for i in $(ls */*txt); do
	echo $AT $i
	grep $AT $i > /dev/null
	[[ $? == 0 ]] && {
	    grep -v ${AT} ${i} > tmp ;
	diff tmp $i ;
	mv tmp $i ;
	}
    done
done
