#!/bin/bash

(
    # for isotope_str in 0nu_Ecut 2nu_full_Ecut tl208_Ecut bi214_Ecut radon_Ecut
    # do
    #     echo "Isotope : $isotope_str"

    #     isotope=$isotope_str

    #     sed -i -e 's@.*TString isotope.*@TString isotope = "'$isotope_str'";@g' classification_application.C

    #     root -l -q classification_application.C

    # done

    sed -i -e 's@.*bool counts.*@bool counts = false;@g' bdt_score_cc.C

    root -l -b -q bdt_score_cc.C

    sed -i -e 's@.*bool counts.*@bool counts = true;@g' bdt_score_cc.C

    root -l -b -q bdt_score_cc.C

)
