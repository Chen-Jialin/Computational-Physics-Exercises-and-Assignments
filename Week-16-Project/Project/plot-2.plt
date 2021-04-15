set term epslatex standalone 12
set grid
# unset key
set xlabel '$\ln(T_c-T)$ / a.u.'
set ylabel '$\ln m$ / a.u.'
y(x) = a * x + b
fit y(x) 'data-128-128-Wolff-hybrid.txt' u (log(2.27 - $1)):(log($2)) every ::217::226 via a,b
set output '128-128-fit-m-T.tex'
plot 'data-128-128-Wolff-hybrid.txt' u (log(2.27 - $1)):(log($2)) every ::217::226 w p pt 7 ps 1 t 'data point',\
y(x) w l lw 2 t sprintf('fitting line: $\ln m=$%.3f$\ln(2.27-T)+$%.3f',a,b)
set output
set xlabel '$\ln(T_c-T)$ / a.u.'
set ylabel '$\ln\chi$ / a.u.'
fit y(x) 'data-128-128-Wolff-hybrid.txt' u (log(2.27 - $1)):(log($3/128/128)) every ::215::224 via a,b
set output '128-128-fit-chi-T.tex'
plot 'data-128-128-Wolff-hybrid.txt' u (log(2.27 - $1)):(log($3/128/128)) every ::215::224 w p pt 7 ps 1 t 'data point',\
y(x) w l lw 2 t sprintf('fitting line: $\ln\chi=$%.2f$\ln(2.27-T)+$%.2f',a,b)
set output
