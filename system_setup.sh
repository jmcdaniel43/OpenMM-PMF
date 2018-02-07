#!/bin/bash

root=..
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

        $root/bin/build_electrode.pl $root/setup_data/graph.pdb $j > electrode.tmp
        $root/bin/deleteAtoms.py pore $i electrode.tmp > electrode.pdb
        rm electrode.tmp

        cat electrode.pdb $root/setup_data/electrolyte.pdb > SC_start_${j}_${i}.pdb
        rm electrode.pdb

        # echo "Deleting remarks from packmol"
        # elLines=`wc -l $root/setup_data/electrolyte.pdb | awk '{print $1}'`
        # totalLines=`wc -l SC_start_${j}_${i}.pdb | awk '{print $1}'`
        # startOfRemark=`expr $totalLines - $elLines - 1` # delete up to the END line of the first file
        # ending=`expr $startOfRemark +  6`
        # sed -i -e "${startOfRemark},${ending}d" SC_start_${j}_${i}.pdb


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
