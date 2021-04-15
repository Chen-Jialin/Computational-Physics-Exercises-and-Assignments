# set term postscript eps enhanced color font 'Helvetica, 20' 
# set output 'band_plt.eps'
set term x11

 set ylabel 'Energy (ev)' rotate by 90
set yrange [-3.000: 3.000]
set xtics ( ' '  0.0, " " 0.155,)



 plot "band.dat" w l lc rgb "#0060ad" notitle
