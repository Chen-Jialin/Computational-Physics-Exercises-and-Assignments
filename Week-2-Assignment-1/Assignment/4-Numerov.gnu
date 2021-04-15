set term epslatex standalone 12
unset key
set xlabel '$x$'
set ylabel '$\psi(x)$'
set grid 
set output '4-Numerov-0.tex'
plot '4-wavefunction-Numerov-0.txt' w l lc -1 lw 2
set output
