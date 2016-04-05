#!/bin/bash

set -x

echo "Clean stuff"
rm -f fort.* 

echo "Run ATLAS9 to produce a model atmosphere."
./atlas9.bash kapp00p00.ros p00p00bigvt2.bdf ap00t5777g44377k1odfnew.dat t05000_g+3.5_m01p04.at9.mod 5000 +3.5 abundances_c14_p00p00.txt

echo "Test if model converged."
./ifconv.pl t05000_g+3.5_m01p04.at9.mod.log

echo "Run SYNTHE to produce a spectrum."
./synthe.csh
