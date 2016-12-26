#!/bin/bash

(
    for isotope_str in 0nu 2nu
    do
        echo "Isotope : $isotope_str"

        isotope=$isotope_str
        bb='2nu'

        sed -i -e 's@.*TString isotope.*@TString isotope = "'$isotope_str'";@g' classification_application.C

        root -l -q classification_application.C

    done
        root -l -b bdt_score_simple.C

)
