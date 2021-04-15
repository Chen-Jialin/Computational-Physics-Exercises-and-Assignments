set term epslatex standalone 12
set grid
# unset key
set xlabel '$\ln|T_{\max}-T_c|$ / a.u.'
set ylabel '$-\ln L$ / a.u.'
y(x) = a * x + b
fit y(x) 'summary-L-vary.txt' u (log($2 - (2 / log(1 + sqrt(2))))):(-log($1)) via a,b
set output 'fit-L-T.tex'
plot 'summary-L-vary.txt' u (log($2 - (2 / log(1 + sqrt(2))))):(-log($1)) w p pt 7 ps 1 t 'data point',\
y(x) w l lw 2 t sprintf('fitting line: $-\ln L=$%.3f$\ln|T_{\max}-T_c|+$%.3f',a,b)
set output