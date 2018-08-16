#!/bin/bash

for solv in {"acn","dce","h2o"}
do
    for sys in {"bf4_bmim","bf4_tma_tmea"}
    do
        for ion in {"BF4","TMA","TMEA"}
        do
            for sheets in {1..1}
            do
                for pore in {7,10,14}
                do
                    if  [ $ion = "TMEA" ] || [ $ion = "TMA" ] && [ $sys = "bf4_bmim" ]; then
                        continue
                    fi

                    if  [ $ion = "TMEA" ] || [ $ion = "TMA" ] && [ $pore = 7 ]; then
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
