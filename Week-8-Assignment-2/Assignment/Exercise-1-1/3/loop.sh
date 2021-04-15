#! /bin/bash
BIN='/home/JSWL/bin/vasp_ncl'
rm WAVECAR SUMMARY.stanene
for i in 4.662 4.663 4.664 4.665 4.666 ; do
for j in 0.5092 0.5093 0.5094 0.5095 0.5096 ; do
cat >POSCAR  <<!
stanene
$i
        1.0000000000         0.0000000000         0.0000000000
       -0.5000000000         0.8660254037         0.0000000000
        0.0000000000         0.0000000000        20.0000000000/$i
   Sn
    2
Direct
     0.333333349         0.666666698         0.500000000
     0.666666675         0.333333349         $j
!
echo "a= $i  z= $j" ; mpirun -np 16 $BIN
E=`awk '/F=/ {print $0}' OSZICAR` ; echo $i $j $E >>SUMMARY.stanene
done
echo "" >>SUMMARY.stanene
done
cat SUMMARY.stanene

