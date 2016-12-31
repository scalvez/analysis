#!/bin/bash

(
        make

        ./collie_io_bdt.exe

        ./limit_calc.exe limit_bdt.root io_file_bdt.root

        if grep 'Scaling factor' log_collie_bdt.txt; then
            scaling_factor=$(grep 'Scaling factor' log_collie_bdt.txt | awk '{print $4}')
        fi
        echo " Scaling factor :" && echo $scaling_factor

)
