#!/bin/bash

(
    # for sample in confA_A confA_B confA_C confA_D confB_A confB_B confB_C \
    #               confB_D confC_A confC_B confC_C confC_D confD_A confD_B \
    #               confD_C confD_D
    for sample in confA_A confA_B confA_C confA_D confB_A confB_B confB_C \
                  confB_D confC_A confC_B confC_C confC_D
    do
        echo "sample : $sample"
        sed -i -e \
            's@.*cfile->initFile.*@  cfile->initFile("./collie_io/cc/io_file_bdt_'$sample'.root", "0nu");@g' collie_io_bdt.cc

        sed -i -e \
            's@.*TFile infile.*@  TFile infile("./bdt_scores/cc/0nu_rhc/bdt_scores_counts_'$sample'.root");@g' collie_io_bdt.cc

        make
        ./collie_io_bdt.exe
    done
)
