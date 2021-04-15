# set term postscript eps enhanced color font 'Helvetica, 20' 
# set output 'band_plt.eps'
# set term x11

set term epslatex standalone 12

set ylabel 'Energy (ev)' rotate by 90
set yrange [-0.600: 0.600]
set xtics ("$X'$" -0.145,'$\Gamma$'  0.0, "$X'$" 0.145)

set arrow nohead from -0.145,0 to 0.145,0 dt 2 lc rgb 'gray'

set output '1-3-stanene-F.tex'

plot "band-symmetric.dat" w l lw 6 lc rgb 'black' notitle

set output
