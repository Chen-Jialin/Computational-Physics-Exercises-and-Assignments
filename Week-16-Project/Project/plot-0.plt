set term epslatex standalone 12
# set grid
# unset key
set xlabel '$x$'
set ylabel '$y$'
set xrange [0:75]
set yrange [0:75]
set output 'lattice-1.tex'
plot 'lattice-mpi-1.txt' u 1:($3==1?$2:1/0) w p pt 4 t 'spin up', '' u 1:($3==-1?$2:1/0) w p pt 5 t 'spin down'
set output
set output 'lattice-2.tex'
plot 'lattice-mpi-2.txt' u 1:($3==1?$2:1/0) w p pt 4 t 'spin up', '' u 1:($3==-1?$2:1/0) w p pt 5 t 'spin down'
set output
set output 'lattice-3.tex'
plot 'lattice-mpi-3.txt' u 1:($3==1?$2:1/0) w p pt 4 t 'spin up', '' u 1:($3==-1?$2:1/0) w p pt 5 t 'spin down'
set output