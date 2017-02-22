#!/bin/bash

(
    # for sample in confA_A confA_B confA_C confA_D confB_A confB_B confB_C \
    #               confB_D confC_A confC_B confC_C confC_D confD_A confD_B \
    #               confD_C confD_D
    # for sample in confA_A confA_B confA_C confA_D confB_A confB_B confB_C \
    #               confB_D confC_A confC_B confC_C confC_D
    for sample in confD_A confD_B confD_C confD_D
    do
        echo "sample : $sample"
        touch limit_$sample.txt
        time ./limit_calc.exe collie_io/cc/limit_bdt_$sample.root \
        collie_io/cc/io_file_bdt_$sample.root &> limit_$sample.txt

    done
)
