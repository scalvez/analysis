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
            's@.*TFile \* f_0nu.*@TFile \* f_0nu = TFile::Open("./bdt_scores/cc/0nu_rhc/individual_channel/0nu_Ecut_'$sample'.root");@g' bdt_score_cc.C
        sed -i -e \
            's@.*TFile \* f_2nu.*@TFile \* f_2nu = TFile::Open("./bdt_scores/cc/0nu_rhc/individual_channel/2nu_full_Ecut_'$sample'.root");@g' bdt_score_cc.C
        sed -i -e \
            's@.*TFile \* f_tl208.*@TFile \* f_tl208 = TFile::Open("./bdt_scores/cc/0nu_rhc/individual_channel/tl208_Ecut_'$sample'.root");@g' bdt_score_cc.C
        sed -i -e \
            's@.*TFile \* f_bi214.*@TFile \* f_bi214 = TFile::Open("./bdt_scores/cc/0nu_rhc/individual_channel/bi214_Ecut_'$sample'.root");@g' bdt_score_cc.C
        sed -i -e \
            's@.*TFile \* f_radon.*@TFile \* f_radon = TFile::Open("./bdt_scores/cc/0nu_rhc/individual_channel/radon_Ecut_'$sample'.root");@g' bdt_score_cc.C

        sed -i -e \
            's@.*TFile \* f_output= new TFile.*bdt_scores_counts.*@TFile \* f_output= new TFile("./bdt_scores/cc/0nu_rhc/bdt_scores_counts_'$sample'.root","RECREATE");@g' bdt_score_cc.C
        sed -i -e \
            's@.*TFile \* f_output= new TFile.*bdt_scores_conf.*@TFile \* f_output= new TFile("./bdt_scores/cc/0nu_rhc/bdt_scores_'$sample'.root","RECREATE");@g' bdt_score_cc.C

        sed -i -e 's@.*bool counts.*@bool counts = false;@g' bdt_score_cc.C
        root -l -b -q bdt_score_cc.C

        sed -i -e 's@.*bool counts.*@bool counts = true;@g' bdt_score_cc.C
        root -l -b -q bdt_score_cc.C

    done
)
