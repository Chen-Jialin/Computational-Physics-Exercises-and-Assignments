# set term postscript eps enhanced color font 'Helvetica, 20' 
# set output 'band_plt.eps'
set term x11

set ylabel 'Energy (ev)' rotate by 90
set yrange [-15.000: 15.000]
set xtics ( 'L'  0.0, "G" 1.395,"X" 3.006,"U/K" 3.576,"G" 6.271)


set arrow nohead from 1.395,-15.000 to 1.395, 15.000 dt 2 lc rgb "black"
set arrow nohead from 3.006,-15.000 to 3.006, 15.000 dt 2 lc rgb "black"
set arrow nohead from 3.576,-15.000 to 3.576, 15.000 dt 2 lc rgb "black"

 plot "band.dat" w l lc rgb "#0060ad" notitle
