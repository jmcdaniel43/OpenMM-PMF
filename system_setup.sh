#!/bin/bash

root=..
sys="bf4_tmea_dce"

mkdir -p "ffdir"
mkdir -p "pdb"

for j in {1,2}
do
    for i in {7,10,14}
    do
        cd "ffdir"
        $root/bin/deleteAtoms.py pore $i $root/setup_data/graph_name.xml > graph_${i}.xml
        $root/bin/deleteAtoms.py pore $i $root/setup_data/graph_residue.xml > graph_residue_${i}.xml
        cd ..

        cd "pdb"
        cd $sys

        $root/../bin/build_electrode.pl $root/../setup_data/graph.pdb $j > electrode.tmp
        $root/../bin/deleteAtoms.py pore $i electrode.tmp > electrode.pdb
        rm electrode.tmp

        cat electrode.pdb $root/../setup_data/${sys}.pdb > SC_start_${j}_${i}.pdb
        rm electrode.pdb

        cd ..
        cd ..
    done

done

# mkdir -p slurm_scripts
# cd slurm_scripts
# for i in {1..2}
# do
#     for j in {7,10,14}
#     do
#         sed -e "s/XXX/${i}/g" ../run.pbs > run${i}_${j}.pbs
#         sed -i "s/YYY/${j}/g" run${i}_${j}.pbs
#     done
# done
# cd ..
