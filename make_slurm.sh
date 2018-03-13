#!/bin/bash

mkdir -p slurm_scripts
cd slurm_scripts
for solv in {"acn","dce"}
do
    for sys in {"bf4_tmea","bf4_tma_tmea"}
    do
        for ion in {"bf4","tma","tmea"}
        do
            for sheets in {1..2}
            do
                for pore in {7,10,14}
                do
                    scriptName="${sys}_${solv}_${sheets}_${pore}_${ion}diff.pbs"
                    echo $scriptName
                    sed -e "s/SOLVENT/${solv}/g" ../run.pbs > "$scriptName"
                    sed -i "s/SYSTEM/${sys}/g" $scriptName
                    sed -i "s/SHEETS/${sheets}/g" $scriptName
                    sed -i "s/PORE/${pore}/g"   $scriptName
                    sed -i "s/ION/${ion}/g" $scriptName 
                done
            done
        done
    done
done
cd ..
