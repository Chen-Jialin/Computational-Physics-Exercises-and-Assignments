set term epslatex standalone 12
unset key
set xlabel '$x$'
set ylabel '$\Psi(x)$'
set grid
set output '1-2-1.tex'
plot '1-wavefunction.txt' w l lc -1 lw 2
set output
