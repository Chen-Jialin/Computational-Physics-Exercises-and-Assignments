set term epslatex standalone 12
unset key
set xrange [0:1000000]
set xtics ('$0$' 0, '$2\times10^5$' 200000, '$4\times10^5$' 400000, '$6\times10^5$' 600000, '$8\times10^5$' 800000, '$10^6$' 1000000)
set xlabel '$n$'
set ylabel '$F=\int_0^3e^x\,dx=\frac{b-a}{n}\sum_{i=1}^nf(x_i)$'
set arrow nohead from 0,19.08553692 to 1000000,19.08553692 dt 2 lw 4 lc rgb 'black'
set output '1-I-n.tex'
plot '1-I-n.txt'
set output