#!/bin/bash
# this script reads parameter from command line
# run with bash -x atlas9.sh for debugging

# $1 rosseland file
# $2 odf file
# $3 initial model
# $4 output file
# $5 teff
# $6 logg
# $7 abundance file

set -e

echo "============= Starting computation of "$4
/bin/rm -f ./fort.*
if [ -f $4 ]; then #output file exists
   echo 'Output model already exists'
   exit 1
fi

echo 'rosseland file ' $1
echo 'odf ' $2
echo 'initial model' $3
echo 'output ' $4
echo 'teff ' $5
echo 'logg ' $6
echo 'abundances ' $7

export KURUCZ=~/Documents/kurucz_light/
rm -f fort.*
date

abundances="`cat $7`"

#Opacity and ODF file calls
ln -s  $KURUCZ/ODF/Coelho14/$1          fort.1   #kapp00.ros
ln -s  $KURUCZ/ODF/Coelho14/$2          fort.9   #p00big1.bdf
ln -s  $KURUCZ/molecules/molecules.dat fort.2

#Starting model assignation
ln -s $3  fort.3        #t4904g240p00a0vt10_l.mod

#Atlas is called and fed with his input control cards
$KURUCZ/bin/atlas9mem.exe <<EOF>$4.log
READ KAPPA
READ PUNCH
MOLECULES ON
READ MOLECULES
FREQUENCIES 337 1 337 BIG
VTURB 2.0E+5
CONVECTION OVER 1.25 0 36
TITLE ATLAS9 $5 $6 $7
echo $abundances
ITERATIONS 15  PRINT 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
PUNCH 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
SCALE 72 -6.875 0.125 $5 $6
BEGIN                    ITERATION  10 COMPLETED
ITERATIONS 15  PRINT 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
PUNCH 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
SCALE 72 -6.875 0.125 $5 $6
BEGIN                    ITERATION  10 COMPLETED
END
EOF

#The exit model is renamed
mv fort.7 $4

#Some garbage is removed
rm -f fort.1
rm -f fort.9
rm -f fort.2
rm -f fort.3
date
