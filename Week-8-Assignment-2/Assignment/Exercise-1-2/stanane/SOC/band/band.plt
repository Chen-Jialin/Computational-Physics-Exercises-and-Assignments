# set term postscript eps enhanced color font 'Helvetica, 20' 
# set output 'band_plt.eps'
# set term x11

set ylabel 'Energy (ev)' rotate by 90
set yrange [-1.800: 1.800]
set xtics ('M'  0.0, '$\Gamma$' 0.766, '$K$' 1.650, '$M$' 2.092,)

set arrow nohead from 0.766,-1.800 to 0.766, 1.800 lc rgb 'gray'
set arrow nohead from 1.650,-1.800 to 1.650, 1.800 lc rgb 'gray'

set arrow nohead from 0,0 to 2.091641,0 dt 4 lc rgb 'black'

plot "band.dat" w l lw 4 lc rgb 'red' t 'SOC'
