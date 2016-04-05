#!/bin/bash

set -x

echo "Cleaning fort* files"
/bin/rm -f fort.* 

#echo "Running ATLAS9 to produce a model atmosphere."
#./atlas9.bash kapp00p00.ros p00p00bigvt2.bdf ap00t5777g44377k1odfnew.dat t08200_g+4.25_p00at9.mod 8200 +4.25 abundances_c14_p00p00.txt

#echo "Test if model converged."
#./ifconv.pl t06650_g+4.0_p00at9.mod.log

cp -f in.sun.noCr in
echo "1. Run SYNTHE to produce solar spectrum."
./synthe_trunc.csh

cp -f in6650.noCr in
echo "2. Run SYNTHE to produce Procyon spectrum."
./synthe_trunc.csh

cp -f in8200.noCr in
echo "3. Run SYNTHE to produce Beta Pictoris spectrum."
./synthe_trunc.csh

