#!/bin/bash

(
    for isotope_str in 0nu 2nu
    do
        echo "Isotope : $isotope_str"

        if [ "$isotope_str" = "0nu" ]; then
            sed -i -e 's@.*TString isotope.*@TString isotope = "'$isotope_str'_1M";@g' classification_application_simple.C
       fi

        if [ "$isotope_str" = "2nu" ]; then
            sed -i -e 's@.*TString isotope.*@TString isotope = "'$isotope_str'_full_1M";@g' classification_application_simple.C
       fi

        root -l -q classification_application_simple.C

    done

    sed -i -e 's@.*bool counts.*@bool counts = false;@g' bdt_score_simple.C

    root -l -b bdt_score_simple.C

    sed -i -e 's@.*bool counts.*@bool counts = true;@g' bdt_score_simple.C

    root -l -b bdt_score_simple.C

)
