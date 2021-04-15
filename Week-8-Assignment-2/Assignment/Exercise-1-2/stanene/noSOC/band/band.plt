# set term postscript eps enhanced color font 'Helvetica, 20' 
# set output 'band_plt.eps'
# set term x11

set ylabel 'Energy (ev)' rotate by 90
set yrange [-1.800: 1.800]
set xtics ('$M$'  0.0, '$\Gamma$' 0.774, '$K$' 1.668, '$M$' 2.115)


set arrow nohead from 0.774,-1.800 to 0.774, 1.800 lc rgb 'gray'
set arrow nohead from 1.668,-1.800 to 1.668, 1.800 lc rgb 'gray'

set arrow nohead from 0,0 to 2.114773,0 dt 2 lc rgb 'black'

plot 'band.dat' w l lw 2 dt 2 lc rgb 'black' t 'noSOC'
