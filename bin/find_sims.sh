#!/bin/bash

for path in "$@"
do
	did=`tail $path | grep "Done!" | wc -l`
	if [ "$did" -gt 0 ]
	then
		echo $path | sed -e "s/output_logs\///g"
	fi

done
