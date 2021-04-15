# set term postscript eps enhanced color font 'Helvetica, 20' 
# set output 'band_plt.eps'
# set term x11

set term epslatex standalone 12

set ylabel 'Energy (ev)' rotate by 90
set yrange [-0.600: 1.000]
set xtics ("$-X'$" -0.153, '$\Gamma$'  0.0, "$X'$" 0.153)

set arrow nohead from 0,-0.6 to 0,1 lc rgb 'gray'
set arrow nohead from -0.153,0 to 0.153,0 dt 2 lc rgb 'gray'

set output '1-3-stanane.tex'

plot 'band-symmetric.dat' w l lw 4 lc rgb 'black' notitle

set output
