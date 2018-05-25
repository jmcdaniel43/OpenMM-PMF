#!/bin/bash

pmf_files=`find $1 | grep --color=none ".*diff/pmf$"`
rdf_files=`find $1 | grep --color=none ".*diff/rdf.log$"`
outpath=./$2


for i in $pmf_files
do
    bin/make_grace.py pmf $i/pmf > ${outpath}/pmf_`echo $i | bin/simpath2outputname.py`
done

for i in $rdf_files
do
    bin/make_grace.py rdf $i > ${outpath}/rdf_`echo $i | bin/simpath2outputname.py`
done
