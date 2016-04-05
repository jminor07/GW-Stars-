#!/bin/tcsh
#output
#model atmosphere
#abundances
#tio/notio

set ATLAS_DIR=~/Documents/kurucz_light/

#-----------------------------------------------------------------------------------------

set p1=`sed -n -e 1p in`
set p2=`sed -n -e 2p in`
set p3=`sed -n -e 3p in`
set p4=`sed -n -e 4p in`

echo "i read the following parameters: "
echo "output: "$p1
echo "model atmosphere: "$p2
echo "abundances: "$p3
echo "tio? "$p4
#-----------------------------------------------------------------------------------------
date
echo "initializing..." 
rm -f fort.*
cat synthe.header > synthe.mod
head -n 3 $p2 >> synthe.mod
cat conv.txt >> synthe.mod 
cat $p3 >> synthe.mod
set num = `cat $p2 | wc -l`
echo "number of lines in model atmosfere: $num"
if ($num == 97) then                #atlas9
   tail -n 75 $p2 >> synthe.mod
endif
if ($num == 156) then                #marcs
   tail -n 59 $p2 >> synthe.mod
endif
if ($num == 119) then                #atlas12
   tail -n 75 $p2 >> synthe.mod
endif
set model=synthe.mod 
set outspec=$p1
set brspec=br_${outspec}
ln -s $ATLAS_DIR/lines/he1tables.dat fort.18
ln -s $ATLAS_DIR/molecules/molecules.dat fort.2
ln -s $ATLAS_DIR/lines/continua.dat fort.17
#-----------------------------------------------------------------------------------------
date
echo "xnfpelsyn..." 
$ATLAS_DIR/bin/xnfpelsyn.exe < ${model} > /dev/null
#-----------------------------------------------------------------------------------------
date
echo "synbeg..."
$ATLAS_DIR/bin/synbeg.exe <<EOF > synbeg.out
AIR         249.5     700.5    300000.     0.00    0    -10    .001      0    00
AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF        NREAD
EOF
#-----------------------------------------------------------------------------------------
date
echo "rgfalllinesnew..."
ln -s $ATLAS_DIR/lines/gfall.coelho14.dat fort.11
$ATLAS_DIR/bin/rline2.exe > /dev/null
rm -f fort.11
#-----------------------------------------------------------------------------------------
date
echo "rmolecasc..."
ln -s $ATLAS_DIR/molecules/c2ax.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/c2ba.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/c2da.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/c2ea.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/chax.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/chbx.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/chcx.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/cnax.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/cnbx.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/coax.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/coxx.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/h2bx.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/h2cx.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/mghax.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/mghbx.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/nhax.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/nhca.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/ohnew.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/sihnew.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/sioax.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/sioex.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
ln -s $ATLAS_DIR/molecules/sioxx.dat fort.11
$ATLAS_DIR/bin/rmolecasc.exe > /dev/null
rm -f fort.11
 #-----------------------------------------------------------------------------------------
if ( "$p4" == "tio" ) then
date
echo "tio..."
ln -s $ATLAS_DIR/molecules/schwenke.bin fort.11
ln -s $ATLAS_DIR/molecules/eschwenke.bin fort.48
$ATLAS_DIR/bin/rschwenk.exe >  /dev/null
rm -f fort.11
rm -f fort.48
endif
#-----------------------------------------------------------------------------------------
date
echo "synthe..."
$ATLAS_DIR/bin/synthe_pc.exe >  /dev/null
ln -s ${model} fort.5
cat <<EOF >fort.25
0.0       0.        1.        0.        0.        0.        0.        0.
0.
RHOXJ     R1        R101      PH1       PC1       PSI1      PRDDOP    PRDPOW
EOF
#-----------------------------------------------------------------------------------------
date
echo "spectrv..."
$ATLAS_DIR/bin/spectrv_pc.exe >  /dev/null
mv -f fort.7 ${outspec}
ln -s ${outspec} fort.21
#-----------------------------------------------------------------------------------------
date
echo "broaden..."
$ATLAS_DIR/bin/broaden_pc.exe << EOF >  /dev/null
GAUSSIAN    25000.00RESO
EOF
mv -f fort.22 ${brspec}
rm -f fort.2
ln -s ${brspec} fort.1
#-----------------------------------------------------------------------------------------
date
echo "convert to asc..."
$ATLAS_DIR/bin/converfsynnmtoa.exe > /dev/null
mv -f fort.2 ${outspec}
rm -f ${brspec}
#gzip -f ${outspec}
rm -f fort.*
rm -f synthe.mod
rm -f *.out
echo "done."
date
