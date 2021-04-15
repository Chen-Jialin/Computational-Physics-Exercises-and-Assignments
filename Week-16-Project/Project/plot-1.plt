set term epslatex standalone 12
set grid
unset key
set xlabel '$T / K$'
set ylabel '$m$ / a.u.'
set output '128-128-Wolff-m.tex'
plot 'data-128-128-Wolff-200-2000.txt' u 1:2 w lp lw 2 pt 7 ps 1
set output
set xlabel '$T / K$'
set ylabel '$\chi$ / a.u.'
set output '128-128-Wolff-chi.tex'
plot 'data-128-128-Wolff-200-2000.txt' u 1:($3/128/128) w lp lw 2 pt 7 ps 1
set output
set xlabel '$T / K$'
set ylabel '$C$ / a.u.'
set output '128-128-Wolff-C.tex'
plot 'data-128-128-Wolff-200-2000.txt' u 1:4 w lp lw 2 pt 7 ps 1
set output