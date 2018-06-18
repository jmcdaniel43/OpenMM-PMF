#!/bin/bash

pmf_files=`find $1 | grep --color=none ".*diff/pmf$"`
rdf_files=`find $1 | grep --color=none ".*diff/rdf.dat$"`
ion_files=`find $1 | grep --color=none ".*diff/ion_coordination.dat$"`
bulk_solv=`find $1 | grep --color=none ".*diff/bulk_solvation.dat$"`
bulk_ion=`find $1 | grep --color=none ".*diff/bulk_ion.dat$"`
outpath=./$2


for i in $pmf_files
do
    bin/make_grace.py pmf $i/pmf > "${outpath}/pmf_`echo $i | bin/simpath2outputname.py`"
done

for i in $rdf_files
do
    bin/make_grace.py rdf $i > "${outpath}/solv_`echo $i | bin/simpath2outputname.py`"
done

for i in $ion_files
do
    bin/make_grace.py rdf $i ion > "${outpath}/ion_`echo $i | bin/simpath2outputname.py`"
done

for i in $bulk_solv
do
    bin/make_grace.py rdf $i bulk_solv > "${outpath}/bulk_solv_`echo $i | bin/simpath2outputname.py`"
done

for i in $bulk_ion
do
    bin/make_grace.py rdf $i bulk_ion > "${outpath}/bulk_ion_`echo $i | bin/simpath2outputname.py`"
done
