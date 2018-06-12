#!/bin/bash

for solv in {"acn","dce"}
do
    for sys in {"bf4_tmea","bf4_tma_tmea"}
    do
        for ion in {"BF4","TMA","TMEA"}
        do
            for sheets in {1..2}
            do
                for pore in {7,10,14}
                do
                    if [ $sys = "bf4_tmea" ] && [ $ion = "TMA" ]; then
                        continue
                    fi

                    scriptName="${sys}_${solv}_${sheets}_${pore}_${ion}diff"
                    outputPath=${sheets}sheet/${solv}/${pore}pore/output_${sys}_${ion}diff

                    if [ -f symlinks/$scriptName ]; then
                        echo symlinks/$scriptName exists
                    else
                        echo $scriptName
                        ln -s `pwd`/$outputPath symlinks/$scriptName
                    fi

                    if [ -f $outputPath/output.log ]; then
                        echo symlinks/$scriptName exists
                    else
                        ln -s `pwd`/output_logs/${scriptName}.log $outputPath/output.log
                    fi


                done
            done
        done
    done
done
cd ..
