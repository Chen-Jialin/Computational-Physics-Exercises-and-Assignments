set term epslatex standalone 12
unset key
set xlabel '$x$'
set ylabel '$\psi(x)$'
set yrange [-0.12:0.12]
set grid
set output '4-ExactDiagonalization-3.tex'
plot '4-wavefunction-ExactDiagonalization.txt' u 1:4 w l lc -1 lw 2
set output
