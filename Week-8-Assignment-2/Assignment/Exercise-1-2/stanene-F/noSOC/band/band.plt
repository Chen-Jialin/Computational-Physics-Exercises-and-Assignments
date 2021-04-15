# set term postscript eps enhanced color font 'Helvetica, 20' 
# set output 'band_plt.eps'
set term x11

set ylabel 'Energy (ev)' rotate by 90
set yrange [-1.800: 1.800]
set xtics ('$M$'  0.0, '$\Gamma$' 0.722, '$K$' 1.556, '$M$' 1.973)


set arrow nohead from 0.722,-1.800 to 0.722, 1.800 lc rgb 'gray'
set arrow nohead from 1.556,-1.800 to 1.556, 1.800 lc rgb 'gray'

set arrow nohead from 0,0 to 1.972792,0 dt 2 lc rgb 'black'

plot 'band.dat' w l lw 4 dt 4 lc rgb 'black' t 'noSOC'
