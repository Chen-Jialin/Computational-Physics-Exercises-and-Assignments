# set term postscript eps enhanced color font 'Helvetica, 20' 
# set output 'band_plt.eps'
# set term x11

set term epslatex standalone 12

set ylabel 'Energy (ev)' rotate by 90
set yrange [-1.800: 1.800]
set xrange [0: 1.982]
set xtics ('$M$'  0.0, '$\Gamma$' 0.725, '$K$' 1.563, '$M$' 1.982)


set arrow nohead from 0.725,-1.800 to 0.725, 1.800 lc rgb 'gray'
set arrow nohead from 1.563,-1.800 to 1.563, 1.800 lc rgb 'gray'

set arrow nohead from 0,0 to 1.982,0 dt 2 lc rgb 'black'

set output '1-2-stanene-F.tex'

plot 'band-noSOC.dat' w l lw 4 dt 4 lc rgb 'black' t 'noSOC', 'band-SOC.dat' w l lw 4 lc rgb 'red' t 'SOC'

set output