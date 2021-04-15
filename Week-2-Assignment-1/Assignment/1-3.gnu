set term epslatex standalone 12
set key box
set xlabel '$x$'
set ylabel '$\Psi(x)$'
set grid
set output '1-3-1.tex'
plot '1-wavefunction.txt' w p lc -1 lw 2 t 'calculated',\
'1-wavefunction-exact.txt' w l lc -1 lw 2 t 'exact'
set output
