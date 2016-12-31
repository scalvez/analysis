#!/bin/bash

(
        # make

        # ./collie_io_spectra.exe

        # ./limit_calc.exe limit_spectra.root io_file_spectra.root > log_collie_spectra.txt 2>&1

        if grep 'Scaling factor' log_collie_spectra.txt; then
            scaling_factor=$(grep 'Scaling factor' log_collie_spectra.txt | awk '{print $4}')
        fi
        echo " Scaling factor : $scaling_factor"
)
