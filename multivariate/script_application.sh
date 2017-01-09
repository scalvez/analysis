#!/bin/bash

(
    for isotope_str in 0nu_1M 2nu_full_2MeV tl208 bi214 radon
    do
        echo "Isotope : $isotope_str"

        isotope=$isotope_str
        bb='2nu'

        sed -i -e 's@.*TString isotope.*@TString isotope = "'$isotope_str'";@g' classification_application.C

        root -l -q classification_application.C

    done

    sed -i -e 's@.*bool counts.*@bool counts = false;@g' bdt_score.C

    root -l -b bdt_score.C

    sed -i -e 's@.*bool counts.*@bool counts = true;@g' bdt_score.C

    root -l -b bdt_score.C

)
