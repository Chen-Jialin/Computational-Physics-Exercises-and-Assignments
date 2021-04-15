set term epslatex standalone color 11
set yrange [0.16:0.2]
set view 60,210
set pm3d at b
set xlabel 'a'
set ylabel 'z'
set zlabel 'System energy (eV)' rotate by 90
set output '1-1-1.tex'
splot 'SUMMARY-2.stanene' u 1:2:7 w l notitle
set output
