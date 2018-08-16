#!/bin/bash

root=..
solvent=$1 # "acn"
sys=$2 # "bf4_tma_tmea"

mkdir -p "ffdir"
mkdir -p "pdb"

for j in {1,2}
do
    for i in {10,14}
    do
        cd "ffdir"
            $root/bin/deleteAtoms.py pore $i < $root/construct_system/graph_name.xml > graph_${i}.xml
            $root/bin/deleteAtoms.py pore $i < $root/construct_system/graph_residue.xml > graph_residue_${i}.xml
        cd ..

        cd "pdb"
            mkdir -p $solvent
            cd $solvent
                mkdir -p $sys
                cd $sys

                $root/../../bin/build_electrode.pl $root/../../construct_system/graph.pdb $j > electrode.tmp
                $root/../../bin/deleteAtoms.py pore $i < electrode.tmp > electrode.pdb
                rm electrode.tmp

                cat electrode.pdb $root/../../construct_system/${sys}_${solvent}.pdb > SC_start_${j}_${i}.pdb
                rm electrode.pdb

                cd ..
            cd ..
        cd ..
    done

done
