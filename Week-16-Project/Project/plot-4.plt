set term epslatex standalone 12
set grid
# unset key
set xlabel '$\ln L$ / a.u.'
set ylabel '$\ln|M|$ / a.u.'
y(x) = a * x + b
fit y(x) 'summary-L-vary.txt' u (log($1)):(log($3*$1*$1)) via a,b
set output 'fit-M-L.tex'
plot 'summary-L-vary.txt' u (log($1)):(log($3*$1*$1)) w p pt 7 ps 1 t 'data point',\
y(x) w l lw 2 t sprintf('fitting line: $\ln|M|=$%.3f$\ln L+$%.3f',a,b)
set output
set ylabel '$\ln\chi$ / a.u.'
fit y(x) 'summary-L-vary.txt' u (log($1)):(log($4/$1/$1)) via a,b
set output 'fit-chi-L.tex'
plot 'summary-L-vary.txt' u (log($1)):(log($4/$1/$1)) w p pt 7 ps 1 t 'data point',\
y(x) w l lw 2 t sprintf('fitting line: $\ln\chi=$%.3f$\ln L+$%.3f',a,b)
set output