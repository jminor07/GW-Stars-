 set term x11
set size 3,3
set term postscript eps enhanced color "Times-Roman" 65
set output "CAHKFe.eps"
set encoding iso_8859_1
 set xlabel "<Fe> ({/Times-Roman \305})"
 set ylabel "CaHK ({/Times-Roman \305})"
 set title  "CaHK  vs  <Fe>"
#m(x)=  42.7278+  -6.2715*x
m(x)= 45.3915- 6.7937*x
set yrange[10:24]
set xrange[1.5:2.9]
set nokey


f(x)=  54.4867 - 15.9323*x


#plot "indexm.txt"  index 0 using 14:40:($14-$93):($14+$93):($40-$119):($40+$119) with xyerrorbars title "Virgo Data" lt 1 lw 7 ps 3,\
#"indexm.txt"  index 1 using 14:40:($14-$93):($14+$93):($40-$119):($40+$119) with xyerrorbars title "Virgo Data" lt 7 lw 7 ps 3,\
#"Coma.plot"  index 0 using 14:40:($14-$93):($14+$93):($40-$119):($40+$119) with xyerrorbars title "Coma Data" lt 3 lw 7 ps 3,\
#"tmeasure.out" using 14:40:($14-$93):($14+$93):($40-$119):($40+$119) with xyerrorbars title "Sloan Data" lt 2 lw 7 ps 3,\
#"temp.plot" using 14:40:($14-$93):($14+$93):($40-$119):($40+$119) with xyerrorbars title "Toloba Data" lt 5 lw 3 ps 1,\
# m(x) lt -1 lw 3,\
#"cormodels.new" index 0:6 using 14:40 with linespoints title"Model Grid(Age)" lt 3 lw 3,\
#"cormodels.sort" index 0:8 using 14:40 with linespoints title"Model Grid([Fe/H])" lt 4 lw 3


plot "indexm.txt"  index 0 using (($15+$16)/2):40:(($15+$16)/2-sqrt(($94*$94+$95*$95))/2):(($15+$16)/2+sqrt(($94*$94+$95*$95))/2):($40-$119):($40+$119) with xyerrorbars title "Virgo Data" lt 1 lw 7 ps 3,\
"indexm.txt"  index 1 using (($15+$16)/2):40:(($15+$16)/2-sqrt(($94*$94+$95*$95))/2):(($15+$16)/2+sqrt(($94*$94+$95*$95))/2):($40-$119):($40+$119) with xyerrorbars title "Virgo Data" lt 7 lw 7 ps 3,\
"Coma.plot"  index 0 using (($15+$16)/2):40:(($15+$16)/2-sqrt(($94*$94+$95*$95))/2):(($15+$16)/2+sqrt(($94*$94+$95*$95))/2):($40-$119):($40+$119) with xyerrorbars title "Coma Data" lt 3 lw 7 ps 3,\
"tmeasure.out" using (($15+$16)/2):40:(($15+$16)/2-sqrt(($94*$94+$95*$95))/2):(($15+$16)/2+sqrt(($94*$94+$95*$95))/2):($40-$119):($40+$119) with xyerrorbars title "Sloan Data" lt 2 lw 7 ps 3,\
"temp.plot" using (($15+$16)/2):40:(($15+$16)/2-sqrt(($94*$94+$95*$95))/2):(($15+$16)/2+sqrt(($94*$94+$95*$95))/2):($40-$119):($40+$119) with xyerrorbars title "Tolaba Data" lt 5 lw 3 ps 1,\
 f(x) lt -1 lw 3,\
"cormodels.new" index 0:6 using(($15+$16)/2):40 with linespoints title"Model Grid(Age)" lt 3 lw 3,\
"cormodels.sort" index 0:8 using(($15+$16)/2):40 with linespoints title"Model Grid([Fe/H])" lt 4 lw 3  



pause -1 "\n\t push 'q' and 'return' to exit... \n"